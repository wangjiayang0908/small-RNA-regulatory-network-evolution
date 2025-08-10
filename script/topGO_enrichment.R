# topGO富集分析脚本
# 使用方法:
# Rscript topGO_enrichment.R [工作目录路径]
# 参数说明:
# - 工作目录路径: 包含gene2go文件和基因列表文件的目录路径
# 示例:
# Rscript topGO_enrichment.R /path/to/your/data/directory  # 启用并行计算
# Rscript topGO_enrichment.R /path/to/your/data/directory no_parallel  # 禁用并行计算

# 加载所需的包
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(topGO)
  library(ggplot2)
  library(openxlsx)
  library(gridExtra)
  library(GO.db)
  library(parallel)
  library(doParallel)
  library(foreach)
})

# 设置日志输出函数
log_message <- function(message) {
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"))
}

# 定义工作目录，使用第一个参数作为工作目录
args <- commandArgs(trailingOnly = TRUE)

# 检查参数数量
if (length(args) == 0) {
  stop("错误: 请提供工作目录路径作为第一个参数")
}

work_dir <- args[1]

# 检查是否启用并行计算（第二个参数，可选）
use_parallel <- TRUE
if (length(args) >= 2) {
  if (args[2] == "no_parallel" || args[2] == "FALSE" || args[2] == "false") {
    use_parallel <- FALSE
    cat("用户选择禁用并行计算\n")
  }
}

# 检查工作目录是否存在
if (!dir.exists(work_dir)) {
  stop("错误: 指定的工作目录不存在: ", work_dir)
}

# 设置并行计算参数
setup_parallel <- function() {
  # 检查用户是否选择禁用并行计算
  if (!use_parallel) {
    log_message("用户选择禁用并行计算，使用单线程模式")
    return(NULL)
  }
  
  # 获取系统总核心数
  total_cores <- detectCores()
  
  # 为系统保留2个核心，使用其余核心进行并行计算
  available_cores <- max(1, total_cores - 2)
  
  # 如果可用核心数大于1，则启用并行计算
  if (available_cores > 1) {
    # 注册并行后端
    cl <- makeCluster(available_cores)
    registerDoParallel(cl)
    
    log_message(paste0("启用并行计算，使用 ", available_cores, " 个核心 (总核心数: ", total_cores, ")"))
    
    return(cl)
  } else {
    log_message("可用核心数不足，使用单线程模式")
    return(NULL)
  }
}

# 停止并行计算
stop_parallel <- function(cluster) {
  if (!is.null(cluster)) {
    stopCluster(cluster)
    log_message("已停止并行计算")
  }
}

# 基因组名称
new_names <- c("110R-hap1", "110R-hap2", "5BB-hap1", "5BB-hap2", 
               "MH-hap1", "MH-hap2", "SM-hap1", "SM-hap2")
genomes <- c("110R-hap1", "110R-hap2", "5BB-hap1", "5BB-hap2", 
             "MH-hap1", "MH-hap2", "SM-hap1", "SM-hap2")

# 创建结果存储目录
results_dir <- paste0(work_dir, "results/")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# topGO分析函数
run_topGO_modified <- function(genome, ontology, target_dir, gene2go_file) {
  # 在并行环境中使用cat而不是log_message
  if (exists("log_message")) {
    log_message(paste0("对基因组 ", genome, " 进行 ", ontology, " 本体富集分析"))
  } else {
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 对基因组 ", genome, " 进行 ", ontology, " 本体富集分析\n"))
  }
  
  # 读取目标基因和所有基因列表
  diff_file <- paste0(target_dir, "/", genome, "_all_target_genes_list")
  all_genes_file <- paste0(target_dir, "/", genome, "_all_genes_list")
  
  tryCatch({
    diff <- read.csv(diff_file, header = FALSE)
    diff <- as.character(diff[,1])
    
    all_genes <- read.csv(all_genes_file, header = FALSE)
    all_genes <- as.character(all_genes[,1])
    
    # 导入gene2go数据
    gene2go <- readMappings(gene2go_file, sep = '\t')
    
    # 检查目标基因是否为空
    if(length(diff) == 0) {
      if (exists("log_message")) {
        log_message(paste0("警告: 基因组 ", genome, " 没有目标基因"))
      } else {
        cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 警告: 基因组 ", genome, " 没有目标基因\n"))
      }
      return(data.frame())
    }
    
    # 创建基因列表因子
    genenames <- names(gene2go)
    genelist <- factor(as.integer(genenames %in% diff))
    names(genelist) <- genenames
    
    # 确保基因列表有至少一个位于GO注释中的基因
    if(sum(genelist == 1) == 0) {
      if (exists("log_message")) {
        log_message(paste0("警告: 基因组 ", genome, " 的基因与GO注释不匹配"))
      } else {
        cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 警告: 基因组 ", genome, " 的基因与GO注释不匹配\n"))
      }
      return(data.frame())
    }
    
    # 检查因子级别
    if(length(levels(genelist)) < 2) {
      if (exists("log_message")) {
        log_message(paste0("警告: 因子级别不足，跳过 ", ontology, " 富集分析"))
      } else {
        cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 警告: 因子级别不足，跳过 ", ontology, " 富集分析\n"))
      }
      return(data.frame())
    }
    
    # 创建topGOdata对象
    GOdata <- new("topGOdata", ontology = ontology, allGenes = genelist,
                  annot = annFUN.gene2GO, gene2GO = gene2go)
    
    # 权重算法+Fisher检验 (与CTDGs脚本保持一致)
    test.stat.weight <- new("weightCount", testStatistic = GOFisherTest, name="Fisher test", sigRatio = "ratio")
    resultWeight <- getSigGroups(GOdata, test.stat.weight)
    
    # 经典KS检验
    test.stat.ks <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
    resultKS <- getSigGroups(GOdata, test.stat.ks)
    
    # 运行elim+ks算法
    elim.ks <- runTest(GOdata, algorithm = "elim", statistic = "ks")
    
    # 保留前15个条目
    allRes <- GenTable(GOdata, classic=elim.ks, KS=resultKS, weight = resultWeight,
                       orderBy = "weight", ranksOf = "classic", topNodes = 15)
    
    # 确保结果数据框有足够的行
    if(nrow(allRes) == 0) {
      if (exists("log_message")) {
        log_message(paste0("警告: 基因组 ", genome, " 的 ", ontology, " 富集分析没有找到显著的GO条目"))
      } else {
        cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 警告: 基因组 ", genome, " 的 ", ontology, " 富集分析没有找到显著的GO条目\n"))
      }
      return(data.frame())
    }
    
    # 添加基因组信息和本体类型
    allRes$Group <- genome
    allRes$ontology <- ontology
    
    # 将数值列转换为数值类型
    numeric_cols <- c("Annotated", "Significant", "Expected", "classic", "KS", "weight")
    for(col in numeric_cols) {
      if(col %in% colnames(allRes)) {
        allRes[[col]] <- as.numeric(as.character(allRes[[col]]))
      }
    }
    
    # 获取p值并替换原始统计量
    go_ids <- allRes$GO.ID
    
    # 运行额外的测试以获取p值
    classic_fisher_pvals <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    classic_ks_pvals <- runTest(GOdata, algorithm = "classic", statistic = "ks")
    
    # 提取匹配的p值
    classic_fisher_scores <- score(classic_fisher_pvals)[go_ids]
    classic_ks_scores <- score(classic_ks_pvals)[go_ids]
    
    # 替换到allRes数据框中
    allRes$classic_pvalue <- classic_fisher_scores
    allRes$KS_pvalue <- classic_ks_scores
    
    # 替换原列 - 将classic和KS列的值替换为p值
    allRes$classic <- classic_fisher_scores
    allRes$KS <- classic_ks_scores
    
    return(as.data.frame(allRes))
  }, error = function(e) {
    if (exists("log_message")) {
      log_message(paste0("错误: 在进行基因组 ", genome, " 的 ", ontology, " 富集分析时发生错误: ", e$message))
    } else {
      cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 错误: 在进行基因组 ", genome, " 的 ", ontology, " 富集分析时发生错误: ", e$message, "\n"))
    }
    return(data.frame())
  })
}

# 处理weight值中的范围值并计算对数转换
process_weight_values <- function(data) {
  # 将weight值转换为字符串以便处理
  data$weight <- as.character(data$weight)
  
  # 处理含有"<"符号的范围值
  data$weight <- sapply(data$weight, function(w) {
    if(grepl("<", w)) {
      return(gsub("<", "", w))
    } else {
      return(w)
    }
  })
  
  # 将处理后的weight值转换为数值类型
  data$weight <- as.numeric(data$weight)
  
  # 如果有NA值，替换为一个非常小的默认值
  data$weight[is.na(data$weight)] <- 1e-30
  
  # 计算-log10(weight)值用于点图
  data$weight_log <- -log10(data$weight)
  
  # 对classic和KS列也进行相同的处理
  if("classic" %in% colnames(data)) {
    data$classic <- as.character(data$classic)
    data$classic <- sapply(data$classic, function(c) {
      if(grepl("<", c)) {
        return(gsub("<", "", c))
      } else {
        return(c)
      }
    })
    data$classic <- as.numeric(data$classic)
    data$classic[is.na(data$classic)] <- 1e-30
    data$classic_log <- -log10(data$classic)
  }
  
  if("KS" %in% colnames(data)) {
    data$KS <- as.character(data$KS)
    data$KS <- sapply(data$KS, function(k) {
      if(grepl("<", k)) {
        return(gsub("<", "", k))
      } else {
        return(k)
      }
    })
    data$KS <- as.numeric(data$KS)
    data$KS[is.na(data$KS)] <- 1e-30
    data$KS_log <- -log10(data$KS)
  }
  
  return(data)
}

# 创建GO富集点图
create_go_dotplot <- function(enrichment_data, ont, output_file) {
  log_message(paste0("为 ", ont, " 本体创建点图"))
  
  # 检查数据是否为空
  if(nrow(enrichment_data) == 0) {
    log_message(paste0("警告: 没有 ", ont, " 本体的数据，跳过创建点图"))
    # 创建一个空的点图
    p <- ggplot() + 
      theme_void() + 
      annotate("text", x = 0, y = 0, label = paste0("No significant ", ont, " terms found")) +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    ggsave(output_file, plot = p, width = 8.26, height = 9, units = "in")
    return(p)
  }
  
  # 筛选指定本体的数据
  ont_data <- enrichment_data %>% 
    filter(ontology == ont)
  
  # 确保数值列是数值类型
  numeric_cols <- c("Annotated", "Significant", "Expected", "classic", "KS", "weight")
  for(col in numeric_cols) {
    if(col %in% colnames(ont_data)) {
      ont_data[[col]] <- as.numeric(as.character(ont_data[[col]]))
    }
  }
  
  # 如果筛选后数据为空，则返回一个空的点图
  if(nrow(ont_data) == 0) {
    log_message(paste0("警告: 筛选后没有 ", ont, " 本体的数据，跳过创建点图"))
    p <- ggplot() + 
      theme_void() + 
      annotate("text", x = 0, y = 0, label = paste0("No significant ", ont, " terms found")) +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    ggsave(output_file, plot = p, width = 8.26, height = 9, units = "in")
    return(p)
  }
  
  # 为每个组选择前15个条目
  groups <- unique(ont_data$Group)
  top_results <- data.frame()
  
  for(group in groups) {
    group_data <- ont_data %>% filter(Group == group)
    
    if(nrow(group_data) > 0) {
      # 首先过滤掉Significant值为0的条目
      group_data <- group_data %>% filter(Significant > 0)
      
      if(nrow(group_data) > 0) {
        # 处理weight值中的范围值
        group_data$weight <- as.character(group_data$weight)
        
        # 处理含有"<"符号的范围值
        group_data$weight <- sapply(group_data$weight, function(w) {
          if(grepl("<", w)) {
            return(gsub("<", "", w))
          } else {
            return(w)
          }
        })
        
        # 将处理后的weight值转换为数值类型
        group_data$weight <- as.numeric(group_data$weight)
        
        # 如果有NA值，替换为一个非常小的默认值
        group_data$weight[is.na(group_data$weight)] <- 1e-30
        
        # 按weight值排序并选择前15个条目
        group_top <- group_data %>% 
          arrange(weight) %>% 
          head(15)
        
        # 添加到结果数据框
        top_results <- rbind(top_results, group_top)
      } else {
        log_message(paste0("警告: 基因组 ", group, " 过滤Significant=0后没有数据"))
      }
    }
  }
  
  # 如果没有结果，则返回空的点图
  if(nrow(top_results) == 0) {
    log_message(paste0("警告: 没有选出任何 ", ont, " 本体的条目，跳过创建点图"))
    p <- ggplot() + 
      theme_void() + 
      annotate("text", x = 0, y = 0, label = paste0("No significant ", ont, " terms found")) +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    ggsave(output_file, plot = p, width = 8.26, height = 9, units = "in")
    return(p)
  }
  
  # 处理weight值并计算对数转换
  top_results <- process_weight_values(top_results)
  
  # 将Term转换为因子以保持顺序
  top_results$Term <- factor(top_results$Term, levels = unique(top_results$Term))
  
  # 确保Group列的顺序
  desired_order <- genomes
  valid_groups <- desired_order[desired_order %in% unique(top_results$Group)]
  top_results$Group <- factor(top_results$Group, levels = valid_groups)
  
  # 创建点图
  p <- ggplot(top_results, aes(Term, Group)) +
    geom_point(aes(color = weight_log, size = Significant)) +
    theme_bw() +
    scale_size_continuous(name = "Gene number")+
    scale_color_gradient(low = '#6699CC', high = '#CC3333', name = "-log10(weight)") +
    labs(x = NULL, y = NULL, title = paste0(ont, " Ontology Enrichment")) +
    guides(
      size = guide_legend(order = 1, title.position = "top"),
      color = guide_colorbar(order = 2, title.position = "top")
    ) +
    theme(
      legend.direction = "vertical",
      legend.position = "right",
      legend.title = element_text(size = 7,face = "bold"),
      legend.text = element_text(size = 6),
      legend.key.size = unit(.7, "lines"),
      legend.box = "vertical",
      legend.spacing.y = unit(5, "pt"),
      legend.margin = margin(0, 5, 0, 0),
      axis.text.x = element_text(hjust = 0, size = 7, colour = "black", angle=45),
      axis.text.y = element_text(size = 7, angle=0, color = "black", face = "bold"),
      plot.title = element_blank()
    ) +
    scale_y_discrete(position = "right",limits = rev)+
    scale_x_discrete(position = "top",limits = rev)
  
  # 保存点图
  ggsave(output_file, plot = p, width = 8.26, height = 3.8, units = "in", dpi = 300)
  log_message(paste0("点图已保存到: ", output_file))
  
  return(p)
}

# 保存结果到Excel文件（按本体类型分工作表）
save_to_excel <- function(enrichment_data, output_file) {
  log_message(paste0("保存结果到Excel文件: ", output_file))
  
  # 检查数据是否为空
  if(nrow(enrichment_data) == 0) {
    log_message("警告: 没有数据可以保存到Excel文件")
    wb <- createWorkbook()
    addWorksheet(wb, "Info")
    writeData(wb, "Info", data.frame(Message = "No significant GO terms found"))
    saveWorkbook(wb, output_file, overwrite = TRUE)
    return(invisible(NULL))
  }
  
  # 确保数值列是数值类型
  numeric_cols <- c("Annotated", "Significant", "Expected", "classic", "KS", "weight")
  for(col in numeric_cols) {
    if(col %in% colnames(enrichment_data)) {
      enrichment_data[[col]] <- as.numeric(as.character(enrichment_data[[col]]))
    }
  }
  
  # 创建工作簿
  wb <- createWorkbook()
  
  # 为每个本体创建工作表
  ontologies <- c("MF", "BP", "CC")
  
  for(ont in ontologies) {
    # 筛选本体数据
    ont_data <- enrichment_data %>% 
      filter(ontology == ont)
    
    if(nrow(ont_data) > 0) {
      # 添加工作表
      addWorksheet(wb, ont)
      
      # 为每个基因组选择前15个条目
      groups <- unique(ont_data$Group)
      top_results <- data.frame()
      
      for(group in groups) {
        group_data <- ont_data %>% filter(Group == group)
        
        # 筛选出Significant值大于0的条目
        group_data <- group_data %>% filter(Significant > 0)
        
        # 如果过滤后还有数据，则处理weight值和排序
        if(nrow(group_data) > 0) {
          # 处理weight值中的范围值
          group_data$weight <- as.character(group_data$weight)
          
          # 处理含有"<"符号的范围值
          group_data$weight <- sapply(group_data$weight, function(w) {
            if(grepl("<", w)) {
              return(gsub("<", "", w))
            } else {
              return(w)
            }
          })
          
          # 将处理后的weight值转换为数值类型
          group_data$weight <- as.numeric(group_data$weight)
          
          # 如果有NA值，替换为一个非常小的默认值
          group_data$weight[is.na(group_data$weight)] <- 1e-30
          
          # 按weight值排序并选择前15个条目
          group_top <- group_data %>% 
            arrange(weight) %>% 
            head(15)
          
          # 移除额外生成的对数列和p值列
          extra_cols <- c("weight_log", "classic_log", "KS_log", "classic_pvalue", "KS_pvalue") 
          if(any(extra_cols %in% colnames(group_top))) {
            group_top <- group_top[, !names(group_top) %in% extra_cols]
          }
          
          top_results <- rbind(top_results, group_top)
        } else {
          log_message(paste0("警告: 基因组 ", group, " 过滤Significant=0后没有数据"))
        }
      }
      
      # 确保Group列的顺序
      if(nrow(top_results) > 0) {
        desired_order <- genomes
        valid_groups <- desired_order[desired_order %in% unique(top_results$Group)]
        top_results$Group <- factor(top_results$Group, levels = valid_groups)
        top_results$weight <- as.numeric(as.character(top_results$weight))
        top_results <- top_results[order(top_results$Group, top_results$weight), ]
        
        # 写入数据
        writeData(wb, ont, top_results)
        
        # 设置列宽
        setColWidths(wb, ont, cols = 1:ncol(top_results), widths = "auto")
      } else {
        # 添加空的工作表
        writeData(wb, ont, data.frame(Message = paste0("No significant ", ont, " terms found")))
      }
    } else {
      # 添加空的工作表
      addWorksheet(wb, ont)
      writeData(wb, ont, data.frame(Message = paste0("No significant ", ont, " terms found")))
    }
  }
  
  # 添加汇总工作表
  addWorksheet(wb, "Summary")
  
  # 过滤掉Significant为0的数据用于统计
  filtered_enrichment_data <- enrichment_data %>% filter(Significant > 0)
  
  # 计算每个基因组每个本体的GO条目数
  if(nrow(filtered_enrichment_data) > 0) {
    summary_data <- filtered_enrichment_data %>%
      group_by(Group, ontology) %>%
      summarise(Term_Count = n(), .groups = 'drop') %>%
      pivot_wider(names_from = ontology, values_from = Term_Count, values_fill = 0)
    
    # 如果没有summary数据，添加提示信息
    if(nrow(summary_data) == 0) {
      writeData(wb, "Summary", data.frame(Message = "No significant GO terms found"))
    } else {
      # 确保所有本体类型都存在，如果不存在则添加0值列
      for(ont in ontologies) {
        if(!(ont %in% colnames(summary_data))) {
          summary_data[[ont]] <- 0
        }
      }
      
      # 写入汇总数据
      writeData(wb, "Summary", summary_data)
      setColWidths(wb, "Summary", cols = 1:ncol(summary_data), widths = "auto")
    }
  } else {
    writeData(wb, "Summary", data.frame(Message = "No significant GO terms found after filtering Significant=0"))
  }
  
  # 添加包含所有过滤后结果的工作表
  addWorksheet(wb, "All_Results")
  
  if(nrow(filtered_enrichment_data) > 0) {
    # 移除额外生成的对数列和p值列
    extra_cols <- c("weight_log", "classic_log", "KS_log", "classic_pvalue", "KS_pvalue")
    if(any(extra_cols %in% colnames(filtered_enrichment_data))) {
      filtered_enrichment_data <- filtered_enrichment_data[, !names(filtered_enrichment_data) %in% extra_cols]
    }
    
    # 写入所有数据
    writeData(wb, "All_Results", filtered_enrichment_data)
    
    # 自动调整列宽
    setColWidths(wb, "All_Results", cols = 1:ncol(filtered_enrichment_data), widths = "auto")
  } else {
    writeData(wb, "All_Results", data.frame(Message = "No significant GO terms after filtering Significant=0"))
  }
  
  # 保存工作簿
  saveWorkbook(wb, output_file, overwrite = TRUE)
  log_message(paste0("结果已保存到: ", output_file))

  # 输出BP类别功能条目-基因组名称的二维矩阵excel表格
  # 仅保留每个基因组最显著富集的前15个BP条目
  log_message("生成BP类别功能条目-基因组二维矩阵...")
  bp_data <- enrichment_data %>% filter(ontology == "BP", Significant > 0)
  matrix_list <- list()
  for (group in genomes) {
    group_data <- bp_data %>% filter(Group == group)
    if (nrow(group_data) > 0) {
      # 处理weight值
      group_data$weight <- as.character(group_data$weight)
      group_data$weight <- sapply(group_data$weight, function(w) {
        if(grepl("<", w)) {
          return(gsub("<", "", w))
        } else {
          return(w)
        }
      })
      group_data$weight <- as.numeric(group_data$weight)
      group_data$weight[is.na(group_data$weight)] <- 1e-30
      # 取前15个
      group_top <- group_data %>% arrange(weight) %>% head(15)
      matrix_list[[group]] <- group_top[, c("GO.ID", "Term")]
    } else {
      matrix_list[[group]] <- data.frame(GO.ID=character(0), Term=character(0))
    }
  }
  # 合并所有基因组的GO.ID-Term，去重
  all_go_terms <- do.call(rbind, matrix_list)
  all_go_terms <- unique(all_go_terms)
  # 构建矩阵
  matrix_df <- data.frame(GO.ID = all_go_terms$GO.ID, Term = all_go_terms$Term, stringsAsFactors = FALSE)
  for (group in genomes) {
    group_gos <- matrix_list[[group]]$GO.ID
    matrix_df[[group]] <- ifelse(matrix_df$GO.ID %in% group_gos, 1, 0)
  }
  # 输出为新Excel文件
  bp_matrix_file <- sub("by_Ontology\\.xlsx$", "BP_Matrix.xlsx", output_file)
  wb2 <- createWorkbook()
  addWorksheet(wb2, "BP_Matrix")
  writeData(wb2, "BP_Matrix", matrix_df)
  setColWidths(wb2, "BP_Matrix", cols = 1:ncol(matrix_df), widths = "auto")
  saveWorkbook(wb2, bp_matrix_file, overwrite = TRUE)
  log_message(paste0("BP功能条目-基因组矩阵已保存到: ", bp_matrix_file))
}

# 主分析函数
main <- function() {
  log_message("开始运行修改后的topGO富集分析脚本")
  
  # 设置并行计算
  cluster <- setup_parallel()
  
  # 设置文件路径
  target_dir <- work_dir
  gene2go_file <- paste0(target_dir, "/gene2go")
  
  # 创建一个空列表来存储所有结果
  all_results <- list()
  combined_results <- data.frame()
  
  # 准备并行分析的数据
  analysis_tasks <- data.frame(
    genome = genomes,
    new_name = new_names,
    stringsAsFactors = FALSE
  )
  
  # 根据是否启用并行计算选择执行方式
  if (!is.null(cluster)) {
    # 并行执行所有基因组的分析
    log_message("开始并行分析所有基因组...")
    
    # 使用foreach进行并行计算
    results_list <- foreach(i = 1:nrow(analysis_tasks), 
                           .packages = c("topGO", "dplyr", "tidyr", "ggplot2", "openxlsx", "gridExtra", "GO.db"),
                           .export = c("run_topGO_modified", "log_message", "target_dir", "gene2go_file", "analysis_tasks", "process_weight_values"),
                           .combine = "c",
                           .errorhandling = "pass") %dopar% {
      
      genome <- analysis_tasks$genome[i]
      new_name <- analysis_tasks$new_name[i]
      
      # 创建日志消息（在并行环境中需要特殊处理）
      cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 分析基因组 ", genome, " (", new_name, ")...\n"))
      
      # 对每个ontology运行分析
      result_BP <- run_topGO_modified(genome, "BP", target_dir, gene2go_file)
      result_MF <- run_topGO_modified(genome, "MF", target_dir, gene2go_file)
      result_CC <- run_topGO_modified(genome, "CC", target_dir, gene2go_file)
      
      # 合并三个ontology的结果
      genome_result <- rbind(result_BP, result_MF, result_CC)
      
      # 返回结果（使用列表包装以便combine）
      list(setNames(list(genome_result), genome))
    }
    
    # 处理并行结果
    for (result_item in results_list) {
      genome_name <- names(result_item)[1]
      genome_result <- result_item[[1]]
      
      # 保存到列表中
      all_results[[genome_name]] <- genome_result
      
      # 添加到组合结果
      combined_results <- rbind(combined_results, genome_result)
    }
  } else {
    # 使用传统的for循环进行串行分析
    log_message("开始串行分析所有基因组...")
    
    for (i in 1:nrow(analysis_tasks)) {
      genome <- analysis_tasks$genome[i]
      new_name <- analysis_tasks$new_name[i]
      
      log_message(paste0("分析基因组 ", genome, " (", new_name, ")..."))
      
      # 对每个ontology运行分析
      result_BP <- run_topGO_modified(genome, "BP", target_dir, gene2go_file)
      result_MF <- run_topGO_modified(genome, "MF", target_dir, gene2go_file)
      result_CC <- run_topGO_modified(genome, "CC", target_dir, gene2go_file)
      
      # 合并三个ontology的结果
      genome_result <- rbind(result_BP, result_MF, result_CC)
      
      # 保存到列表中
      all_results[[genome]] <- genome_result
      
      # 添加到组合结果
      combined_results <- rbind(combined_results, genome_result)
    }
  }
  
  # 检查是否有结果
  if(nrow(combined_results) > 0) {
    # 设置输出文件路径
    output_xlsx <- file.path(results_dir, "GO_Enrichment_by_Ontology.xlsx")
    output_mf_plot <- file.path(results_dir, "GO_Enrichment_MF.pdf")
    output_bp_plot <- file.path(results_dir, "GO_Enrichment_BP.pdf")
    output_cc_plot <- file.path(results_dir, "GO_Enrichment_CC.pdf")
    
    # 创建可视化
    log_message("创建可视化")
    mf_plot <- create_go_dotplot(combined_results, "MF", output_mf_plot)
    bp_plot <- create_go_dotplot(combined_results, "BP", output_bp_plot)
    cc_plot <- create_go_dotplot(combined_results, "CC", output_cc_plot)
    
    # 保存结果到Excel文件
    save_to_excel(combined_results, output_xlsx)
    
    log_message(paste0("处理完成，结果已保存到目录: ", results_dir))
  } else {
    log_message("警告: 没有找到任何显著富集的GO条目")
    
    # 创建空的点图
    for(ont in c("MF", "BP", "CC")) {
      output_file <- file.path(results_dir, paste0("GO_Enrichment_", ont, ".pdf"))
      p <- ggplot() + 
        theme_void() + 
        annotate("text", x = 0, y = 0, label = paste0("No significant ", ont, " terms found")) +
        theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
      ggsave(output_file, plot = p, width = 8.26, height = 9, units = "in")
    }
    
    # 创建空的Excel文件
    output_xlsx <- file.path(results_dir, "GO_Enrichment_by_Ontology.xlsx")
    wb <- createWorkbook()
    addWorksheet(wb, "Info")
    writeData(wb, "Info", data.frame(Message = "No significant GO terms found"))
    saveWorkbook(wb, output_xlsx, overwrite = TRUE)
    
    log_message(paste0("处理完成，但没有找到显著结果，已创建空白输出文件: ", results_dir))
  }
  
  log_message("脚本执行完成")
  
  # 停止并行计算
  stop_parallel(cluster)
}

# 执行主函数
main()