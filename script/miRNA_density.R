#!/usr/bin/env Rscript
# miRNA前体染色体密度分析脚本
# 使用滑动窗口分析8个基因组中miRNA前体的密度分布

# 加载必要的R包
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(viridis)
library(RColorBrewer)
library(stringr)

setwd("~/Desktop/project/sRNA_prediction/比较基因组学解析鲜食和砧木葡萄群体间 miRNA 靶基因/data for visualization/sliding_window_data")
#################################################
# 参数设置
#################################################

# 基因组列表
genomes <- c("110R-hap1", "110R-hap2", "5BB-hap1", "5BB-hap2", "MH-hap1", "MH-hap2", "SM-hap1", "SM-hap2")

# 文件路径模板
gff_template <- "%s_miRNA_precursor.gff"  # miRNA前体GFF文件
chr_template <- "%s.chr"                   # 染色体长度文件

# 滑动窗口参数
window_size <- 1000000  # 窗口大小 (1Mb)
step_size <- 200000     # 滑动步长 (200kb)

# 输出文件设置
output_dir <- "density_plots"
output_file <- "miRNA_precursor_density.pdf"
output_width <- 8.3
output_height <- 10
output_dpi <- 300

# 颜色设置 - 为8个基因组分配不同颜色
genome_colors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#cab2d6")
names(genome_colors) <- genomes

#################################################
# 函数设置
#################################################

# 读取染色体长度文件
read_chromosome_lengths <- function(genome) {
  file_path <- sprintf(chr_template, genome)
  
  # 检查文件是否存在
  if (!file.exists(file_path)) {
    warning(paste("染色体长度文件不存在:", file_path))
    return(NULL)
  }
  
  # 读取染色体长度文件
  chr_data <- read.table(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(chr_data) <- c("chromosome", "length")
  
  # 添加基因组信息
  chr_data$genome <- genome
  
  return(chr_data)
}

# 读取GFF文件并提取miRNA前体信息
read_mirna_gff <- function(genome) {
  file_path <- sprintf(gff_template, genome)
  
  # 检查文件是否存在
  if (!file.exists(file_path)) {
    warning(paste("GFF文件不存在:", file_path))
    return(NULL)
  }
  
  # 读取GFF文件
  gff_data <- read.table(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  # 提取必要的列: 染色体、起始位置、终止位置、ID
  mirna_data <- gff_data %>%
    select(V1, V4, V5, V9) %>%
    rename(chromosome = V1, start = V4, end = V5, attributes = V9)
  
  # 从attributes列中提取ID
  mirna_data$id <- sub(".*ID=([^;]+).*", "\\1", mirna_data$attributes)
  
  # 添加基因组信息
  mirna_data$genome <- genome
  
  # 移除重复条目 (基于ID和位置) - 确保每个miRNA前体只计算一次
  mirna_data <- mirna_data %>%
    distinct(chromosome, start, end, id, .keep_all = TRUE)
  
  return(mirna_data)
}

# 使用滑动窗口计算miRNA前体密度
calculate_density <- function(mirna_data, chr_lengths, window_size, step_size) {
  # 初始化结果数据框
  density_results <- data.frame()
  
  # 确保染色体排序一致
  all_chromosomes <- unique(chr_lengths$chromosome)
  chr_nums <- as.numeric(gsub("chr([0-9]+).*", "\\1", all_chromosomes))
  if (!all(is.na(chr_nums))) {
    all_chromosomes <- all_chromosomes[order(chr_nums)]
  }
  
  # 为每个基因组和染色体组合计算密度
  for (genome in unique(chr_lengths$genome)) {
    for (chr in all_chromosomes) {
      # 检查当前染色体是否存在于此基因组
      if (!chr %in% chr_lengths$chromosome[chr_lengths$genome == genome]) {
        next
      }
      
      # 获取当前染色体长度
      chr_length <- chr_lengths$length[chr_lengths$genome == genome & 
                                         chr_lengths$chromosome == chr]
      
      # 过滤当前基因组和染色体的miRNA数据
      current_data <- mirna_data %>%
        filter(genome == !!genome, chromosome == !!chr)
      
      # 如果没有miRNA数据，创建一条零密度记录并继续
      if (nrow(current_data) == 0) {
        density_results <- rbind(density_results, data.frame(
          genome = genome,
          chromosome = chr,
          window_start = 1,
          window_end = window_size,
          window_center = window_size/2,
          count = 0,
          density = 0
        ))
        next
      }
      
      # 定义窗口起始位置
      window_starts <- seq(1, max(1, chr_length - window_size + 1), by = step_size)
      if (length(window_starts) == 0) {
        window_starts <- 1
      }
      
      # 对每个窗口计算密度
      for (start_pos in window_starts) {
        end_pos <- min(start_pos + window_size - 1, chr_length)
        actual_window_size <- end_pos - start_pos + 1
        
        # 计算窗口内的miRNA数量
        window_count <- sum(
          current_data$start >= start_pos & current_data$end <= end_pos
        )
        
        # 计算密度 (每百万碱基对/Mb的miRNA数量)
        # 使用实际窗口大小防止边界处窗口不完整
        density <- (window_count / (actual_window_size / 1000000))
        
        # 添加到结果
        density_results <- rbind(density_results, data.frame(
          genome = genome,
          chromosome = chr,
          window_start = start_pos,
          window_end = end_pos,
          window_center = start_pos + actual_window_size/2,
          count = window_count,
          density = density
        ))
      }
    }
  }
  
  return(density_results)
}

###############################################################################
# 创建染色体密度图的函数（曼哈顿图风格，柱状图）
# 染色体密度图函数（按基因组和染色体分网格
###############################################################################
plot_chromosome_density <- function(density_data, genome_colors, chr_lengths, output_file) {
  # 确保包含所有预期的染色体（从chr_lengths中获取）
  all_chromosomes <- unique(chr_lengths$chromosome)
  # 按数字排序染色体
  chr_nums <- as.numeric(gsub("chr([0-9]+).*", "\\1", all_chromosomes))
  if (!all(is.na(chr_nums))) {
    all_chromosomes <- all_chromosomes[order(chr_nums)]
  }
  
  # 获取所有基因组
  all_genomes <- unique(chr_lengths$genome)
  
  # 创建完整的基因组-染色体组合网格
  complete_grid <- expand.grid(
    genome = all_genomes,
    chromosome = all_chromosomes,
    stringsAsFactors = FALSE
  )
  
  # 检查哪些组合在密度数据中缺失
  density_data_keys <- paste(density_data$genome, density_data$chromosome, sep = ":")
  complete_grid_keys <- paste(complete_grid$genome, complete_grid$chromosome, sep = ":")
  missing_combinations <- complete_grid[!complete_grid_keys %in% density_data_keys, ]
  
  if(nrow(missing_combinations) > 0) {
    cat(sprintf("发现 %d 个基因组-染色体组合在密度数据中缺失\n", nrow(missing_combinations)))
    
    # 为每个缺失的组合获取染色体长度
    missing_combinations$window_start <- 1
    missing_combinations$window_end <- window_size
    missing_combinations$window_center <- window_size/2
    missing_combinations$count <- 0
    missing_combinations$density <- 0
    
    # 合并到原始数据
    density_data <- rbind(density_data, missing_combinations)
  }
  
  # 计算最大密度值，确保所有子图有相同的Y轴刻度
  max_density <- max(density_data$density, na.rm = TRUE) * 1.05  # 增加5%空间
  
  # 确保基因组和染色体的正确顺序
  density_data$genome <- factor(density_data$genome, levels = genomes)
  density_data$chromosome <- factor(density_data$chromosome, levels = all_chromosomes)
  
  # 准备绘图数据 - 将染色体位置转换为Mb
  plot_data <- density_data
  plot_data$window_center_mb <- plot_data$window_center / 1000000  # 转换为Mb
  
  # 获取每个染色体的最大长度（以Mb为单位）
  chr_max_lengths <- chr_lengths %>%
    group_by(chromosome) %>%
    summarize(max_length_mb = max(length, na.rm = TRUE) / 1000000)
  
  # 创建刻度计算函数 - 每10Mb一个刻度
  create_mb_breaks <- function(limits) {
    # 确保最大值是10的倍数
    max_mb <- ceiling(limits[2] / 10) * 10
    # 返回从0到最大值的10Mb间隔刻度
    seq(0, max_mb, by = 10)
  }
  
  # 创建输出目录
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  
  # 绘制基因组x染色体网格图
  p <- ggplot(plot_data, aes(x = window_center_mb, y = density, fill = genome)) +
    # 使用柱状图表示密度
    geom_col(width = (step_size / 1000000) * 0.8) +  # 转换柱宽到Mb
    
    # 按基因组和染色体分面，指定drop=FALSE确保显示所有染色体
    facet_grid(genome ~ chromosome, scales = "free_x", space = "free_x", switch = "x", drop = FALSE) +
    
    # 基因组配色
    scale_fill_manual(values = genome_colors) +
    
    # 设置x轴刻度，统一为10Mb间隔
    scale_x_continuous(
      breaks = create_mb_breaks,
      expand = c(0.01, 0)
    ) +
    
    # 固定y轴范围
    scale_y_continuous(limits = c(0, max_density), expand = expansion(mult = c(0, 0))) +
    
    # 设置标题和标签
    labs(
      title = "",
      x = "Position (Mb)",
      y = "Density (miRNAs/Mb)",
      fill = "Genome"
    ) +
    
    # 设置主题
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      panel.spacing = unit(0.2, "lines"),
      strip.text.y = element_text(size = 7, face = "bold", angle = 0),
      strip.text.x = element_text(size = 7, face = "bold", angle = 0),
      strip.placement = "outside",
      strip.background = element_rect(fill = "white", color = "white"),  # 子图背景，移除灰色边框
      axis.text.x = element_text(size = 6, angle = 0,face = "bold"),
      axis.title.x = element_text(size = 12, face="bold"),
      axis.title.y = element_text(size = 12, face="bold"),
      axis.line = element_line(color = "black", size = .5),
      axis.ticks = element_line(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.3), # 子图边框
      plot.margin = margin(1, 1, 1, 1)
    )
  
  # 保存PDF图像
  ggsave(output_file, p, width = output_width, height = output_height, dpi = output_dpi)
  cat(sprintf("密度图已保存为PDF: %s\n", output_file))
  
  # 同时保存PNG图像
  png_output_file <- sub("\\.pdf$", ".png", output_file)
  ggsave(png_output_file, p, width = output_width, height = output_height, dpi = output_dpi)
  cat(sprintf("密度图已保存为PNG: %s\n", png_output_file))
  
  return(p)
}

#################################################
# 主程序
#################################################

# 读取所有基因组的染色体长度信息
chr_length_list <- lapply(genomes, read_chromosome_lengths)
all_chr_lengths <- do.call(rbind, chr_length_list[!sapply(chr_length_list, is.null)])

# 检查是否成功读取染色体长度数据
if (nrow(all_chr_lengths) == 0) {
  stop("未能读取任何染色体长度数据，请检查文件路径和格式")
}

# 读取所有基因组的miRNA前体信息
mirna_data_list <- lapply(genomes, read_mirna_gff)
all_mirna_data <- do.call(rbind, mirna_data_list[!sapply(mirna_data_list, is.null)])

# 检查是否成功读取miRNA数据
if (nrow(all_mirna_data) == 0) {
  stop("未能读取任何miRNA数据，请检查文件路径和格式")
}

# 输出基本统计信息
cat("基本统计信息:\n")
cat(sprintf("共读取 %d 条染色体长度记录\n", nrow(all_chr_lengths)))
cat(sprintf("共读取 %d 条miRNA前体记录\n", nrow(all_mirna_data)))

# 计算密度
density_data <- calculate_density(all_mirna_data, all_chr_lengths, window_size, step_size)

# 输出密度统计信息
cat(sprintf("计算得到 %d 个窗口的密度数据\n", nrow(density_data)))
cat("每个基因组的窗口数量:\n")
print(table(density_data$genome))

#######################
# 创建染色体密度图
######################
output_path <- file.path(output_dir, output_file)
plot <- plot_chromosome_density(density_data, genome_colors, all_chr_lengths, output_path)
