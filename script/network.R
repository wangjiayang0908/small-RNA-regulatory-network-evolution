# miRNA-靶基因网络关系图分析脚本
# 使用igraph包绘制miRNA和其靶基因之间的网络关系

# 加载所需的R包
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(gridExtra)
library(viridis)

# 工作目录
setwd("~/Desktop/project/sRNA_prediction/比较基因组学解析鲜食和砧木葡萄群体间 miRNA 靶基因/data for visualization/network_data/")

#################################################
# 参数设置
#################################################

# 基因组列表
genomes <- c("110R-hap1", "110R-hap2", "5BB-hap1", "5BB-hap2", "MH-hap1", "MH-hap2", "SM-hap1", "SM-hap2")

# 文件路径
file_template <- "%s_all_target_genes_list"

# 布局算法选择 
# 大型网络推荐使用: "drl" 或 "lgl"（用于处理大型网络）
# 其他选项: "fr", "kk", "graphopt", "mds", "sugiyama"
layout_algorithm <- "drl"  # 推荐用于大型网络的布局算法

# 节点大小调整因子
mirna_size_factor <- 5  # miRNA节点大小基础值 (原值: 5)
target_size_factor <- 1 # 靶基因节点大小基础值 (原值: 3)

# 节点大小计算方法
# "log" - 对数缩放 (适合连接度差异大的网络)
# "sqrt" - 平方根缩放 (中等差异)
# "linear" - 线性缩放 (适合连接度相近的网络)
size_scale_method <- "sqrt"

# 最大/最小节点大小限制
max_node_size <- 8    # 节点最大尺寸
min_node_size <- 1    # 节点最小尺寸

# 节点颜色设置
mirna_color <- "#fb8072"  # miRNA节点颜色 (橙色)
target_color <- "#80b1d3"  # 靶基因节点颜色 (蓝色)

# 边线样式
edge_width_val <- 0.3  # 边线宽度 (降低以减少视觉干扰)
edge_color_val <- "#999999"  # 边线颜色
edge_alpha_val <- 0.5  # 边线透明度

# 标签设置
vertex_label_cex_large <- 0.4  # 大型网络中标签字体大小
vertex_label_cex_small <- 0.4  # 小型网络中标签字体大小
hide_labels_threshold <- 500   # 超过此节点数时隐藏标签
show_labels_type <- "all"    # 显示哪类节点的标签: "all", "mirna", "target", "none"

# 单个基因组的简化网络设置
# 为每个基因组创建简化版本以减少节点重叠
genome_top_n_mirna <- 1      # 每个基因组保留的top miRNA数量
genome_top_n_targets <- 100    # 每个基因组保留的top靶基因数量

# 输出图片设置
output_width <- 8.3   # 宽度（英寸）- 增加以容纳更多节点
output_height <- 11.7  # 高度（英寸）
output_dpi <- 400    # 分辨率

#################################################
# 数据读取和处理
#################################################

# 创建一个函数来读取和处理单个基因组的数据
read_mirna_target <- function(genome) {
  file_path <- sprintf(file_template, genome)
  
  # 检查文件是否存在
  if (!file.exists(file_path)) {
    warning(paste("文件不存在:", file_path))
    return(NULL)
  }
  
  # 读取数据
  data <- read.table(file_path, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  colnames(data) <- c("miRNA", "target")
  
  # 为了确保miRNA和靶基因是唯一的，为名称添加基因组前缀
  # 这样可以保证不同基因组间的miRNA和靶基因不会混淆
  data$genome <- genome
  
  return(data)
}

# 读取所有基因组的数据，但保持每个基因组的数据单独存储
genome_data_list <- list()
for (genome in genomes) {
  data <- read_mirna_target(genome)
  if (!is.null(data)) {
    genome_data_list[[genome]] <- data
    cat(paste0("已读取", genome, "的数据，共", nrow(data), "条记录\n"))
  }
}

#################################################
# 创建网络图
#################################################

# 创建网络基本函数
create_network <- function(data, layout_alg="fr") {
  # 创建边列表
  edges <- data.frame(source = data$miRNA, target = data$target)
  
  # 构建图对象
  g <- graph_from_data_frame(edges, directed=TRUE)
  
  # 添加节点类型属性
  V(g)$type <- ifelse(V(g)$name %in% unique(data$miRNA), "miRNA", "target")
  
  # 设置节点大小（基于连接数）- 改进大小缩放方法
  V(g)$degree <- degree(g)
  
  # 计算节点大小，使用选择的缩放方法
  if (size_scale_method == "log") {
    sizes_mirna <- log1p(V(g)$degree[V(g)$type == "miRNA"]) * mirna_size_factor
    sizes_target <- log1p(V(g)$degree[V(g)$type == "target"]) * target_size_factor
  } else if (size_scale_method == "sqrt") {
    sizes_mirna <- sqrt(V(g)$degree[V(g)$type == "miRNA"]) * mirna_size_factor
    sizes_target <- sqrt(V(g)$degree[V(g)$type == "target"]) * target_size_factor
  } else {
    # 线性缩放
    sizes_mirna <- V(g)$degree[V(g)$type == "miRNA"] * mirna_size_factor
    sizes_target <- V(g)$degree[V(g)$type == "target"] * target_size_factor
  }
  
  # 限制节点大小范围
  sizes_mirna <- pmin(pmax(sizes_mirna, min_node_size), max_node_size)
  sizes_target <- pmin(pmax(sizes_target, min_node_size), max_node_size)
  
  # 设置最终的节点大小
  V(g)$size <- 1  # 初始化
  V(g)$size[V(g)$type == "miRNA"] <- sizes_mirna
  V(g)$size[V(g)$type == "target"] <- sizes_target
  
  # 设置节点颜色
  V(g)$color <- ifelse(V(g)$type == "miRNA", mirna_color, target_color)
  
  # 设置节点形状
  V(g)$shape <- ifelse(V(g)$type == "miRNA", "circle", "square")
  
  # 计算布局
  if (layout_alg == "fr") {
    l <- layout_with_fr(g, niter = 500)
  } else if (layout_alg == "kk") {
    l <- layout_with_kk(g)
  } else if (layout_alg == "drl") {
    l <- layout_with_drl(g)
  } else if (layout_alg == "lgl") {
    l <- layout_with_lgl(g)
  } else if (layout_alg == "graphopt") {
    l <- layout_with_graphopt(g, niter = 500)
  } else if (layout_alg == "mds") {
    l <- layout_with_mds(g)
  } else if (layout_alg == "sugiyama") {
    l <- layout_with_sugiyama(g)$layout
  } else {
    l <- layout_with_fr(g, niter = 500)  # 默认使用FR布局
  }
  
  return(list(graph = g, layout = l))
}

#################################################
# 创建基因组特异性的简化网络
#################################################

# 提取高度互连的miRNA和靶基因
create_hub_network <- function(data, top_n_mirna = 10, top_n_targets = 20) {
  # 创建基本网络
  network <- create_network(data, layout_algorithm)
  g <- network$graph
  
  # 找出连接度最高的miRNA节点
  mirna_nodes <- V(g)[V(g)$type == "miRNA"]
  mirna_degrees <- degree(g, v = mirna_nodes)
  
  if (length(mirna_degrees) > 0) {
    top_mirnas <- names(sort(mirna_degrees, decreasing = TRUE)[1:min(top_n_mirna, length(mirna_degrees))])
    
    # 找出与这些miRNA相连的靶基因
    target_nodes <- c()
    for (mirna in top_mirnas) {
      targets <- neighbors(g, mirna, mode = "out")
      target_nodes <- c(target_nodes, targets$name)
    }
    
    # 如果靶基因过多，只保留连接度最高的
    if (length(unique(target_nodes)) > top_n_targets) {
      target_degrees <- degree(g, v = V(g)[V(g)$name %in% unique(target_nodes)])
      top_targets <- names(sort(target_degrees, decreasing = TRUE)[1:top_n_targets])
    } else {
      top_targets <- unique(target_nodes)
    }
    
    # 提取子图
    subgraph_nodes <- c(top_mirnas, top_targets)
    sub_g <- induced_subgraph(g, V(g)[V(g)$name %in% subgraph_nodes])
    
    # 重新计算布局 - 使用FR算法效果更好
    sub_layout <- layout_with_fr(sub_g, niter = 1000)
    
    return(list(graph = sub_g, layout = sub_layout))
  } else {
    return(NULL)
  }
}

# 为每个基因组创建简化网络
hub_networks_by_genome <- list()

for (genome in names(genome_data_list)) {
  cat("为", genome, "创建简化网络...\n")
  hub_networks_by_genome[[genome]] <- create_hub_network(
    genome_data_list[[genome]], 
    top_n_mirna = genome_top_n_mirna, 
    top_n_targets = genome_top_n_targets
  )
}

#################################################
# 绘图函数
#################################################

# 绘制网络图函数 - 支持大型网络的可视化优化
plot_network <- function(network, title = "", vertex_label_cex = 0.6, show_labels = TRUE) {
  g <- network$graph
  layout <- network$layout
  
  # 图例数据
  legend_data <- data.frame(
    type = c("miRNA", "Target Gene"),
    color = c(mirna_color, target_color),
    shape = c(16, 15)  # 16=圆形, 15=方形
  )
  
  # 确定是否显示标签
  display_labels <- show_labels
  if (vcount(g) > hide_labels_threshold) {
    display_labels <- FALSE
    cat("节点数量超过", hide_labels_threshold, "，自动隐藏标签\n")
  }
  
  # 根据设置确定显示哪些标签
  if (display_labels) {
    if (show_labels_type == "mirna") {
      vertex_labels <- ifelse(V(g)$type == "miRNA", V(g)$name, NA)
    } else if (show_labels_type == "target") {
      vertex_labels <- ifelse(V(g)$type == "target", V(g)$name, NA)
    } else if (show_labels_type == "none") {
      vertex_labels <- NA
    } else {
      vertex_labels <- V(g)$name
    }
  } else {
    vertex_labels <- NA
  }
  
  # 绘制基本网络
  plot(g, 
       layout = layout,
       vertex.color = V(g)$color,
       vertex.size = V(g)$size,
       vertex.label = vertex_labels,
       vertex.label.cex = vertex_label_cex,
       vertex.label.dist = 0.5,
       vertex.label.color = "black",
       vertex.shape = ifelse(V(g)$type == "miRNA", "circle", "square"),
       edge.width = edge_width_val,
       edge.color = adjustcolor(edge_color_val, alpha.f = edge_alpha_val),
       edge.arrow.size = 0.2,
       main = title,
       sub = paste("Nodes:", vcount(g), "Edges:", ecount(g), 
                   "| miRNAs:", sum(V(g)$type == "miRNA"), 
                   "Targets:", sum(V(g)$type == "target")),
       frame = FALSE,
       margin = c(0, 0, 0, 0) # 减少边距，增加绘图区域
  )
  
  # 添加图例
  legend("bottomright", 
         legend = legend_data$type,
         col = legend_data$color,
         pch = legend_data$shape,
         pt.cex = 2,
         cex = 0.8,
         bty = "n",
         horiz = FALSE)
}

#################################################
# 保存图片
#################################################

# 创建输出目录
dir.create("network_plots", showWarnings = FALSE)

# 为每个基因组保存简化网络图
for (genome in names(hub_networks_by_genome)) {
  if (!is.null(hub_networks_by_genome[[genome]])) {
    pdf(paste0("network_plots/mirna_target_hub_network_", genome, ".pdf"), 
        width = output_width, height = output_height)
    par(mar = c(1, 1, 3, 1), mfrow = c(1, 1))
    plot_network(hub_networks_by_genome[[genome]], 
                 title = paste("miRNA-Target Gene Hub Network -", genome),
                 vertex_label_cex = vertex_label_cex_small, 
                 show_labels = TRUE)
    dev.off()
    cat("保存", genome, "简化网络图完成\n")
  }
}

# 创建基因组网络比较视图 (每个基因组的简化网络)
if (length(hub_networks_by_genome) > 1) {
  # 过滤掉NULL值
  valid_networks <- hub_networks_by_genome[!sapply(hub_networks_by_genome, is.null)]
  
  if (length(valid_networks) > 0) {
    pdf("network_plots/mirna_target_hub_networks_comparison.pdf", 
        width = output_width * 1.5, height = output_height * 1.5)
    
    # 计算合适的网格大小
    n_plots <- length(valid_networks)
    n_cols <- ceiling(sqrt(n_plots))
    n_rows <- ceiling(n_plots / n_cols)
    
    # 设置网格
    par(mfrow = c(n_rows, n_cols), mar = c(1, 1, 2, 1))
    
    # 绘制每个基因组的简化网络，添加字母标签
    subplot_letters <- letters[1:length(valid_networks)]
    for (i in 1:length(valid_networks)) {
      genome <- names(valid_networks)[i]
      plot_network(valid_networks[[genome]], 
                   title = genome, 
                   vertex_label_cex = 0.4, 
                   show_labels = TRUE)
      
      # 添加字母标签到左上角
      mtext(subplot_letters[i], side = 3, line = -1, adj = 0.05, 
            cex = 1.5, font = 2, col = "black")
    }
    
    dev.off()
    cat("保存基因组简化网络比较图完成\n")
  }
}

#################################################
# 统计信息
#################################################

# 生成单个基因组的统计信息
generate_genome_statistics <- function(data, genome_name) {
  # 计算miRNA数量
  mirna_count <- length(unique(data$miRNA))
  
  # 计算靶基因数量
  target_count <- length(unique(data$target))
  
  # 计算边的数量
  edge_count <- nrow(data)
  
  # 每个miRNA的靶基因数量
  mirna_targets <- table(data$miRNA)
  
  # 计算每个miRNA的平均靶基因数量
  avg_targets <- mean(mirna_targets)
  
  # 找出靶基因数量最多的miRNA
  max_targets_mirna <- names(which.max(mirna_targets))
  max_targets_count <- max(mirna_targets)
  
  # 找出靶基因数量最少的miRNA
  min_targets_mirna <- names(which.min(mirna_targets))
  min_targets_count <- min(mirna_targets)
  
  # 返回统计结果
  return(data.frame(
    Genome = genome_name,
    miRNA_Count = mirna_count,
    Target_Gene_Count = target_count,
    Edge_Count = edge_count,
    Avg_Targets_per_miRNA = round(avg_targets, 2),
    Max_Targets_miRNA = max_targets_mirna,
    Max_Targets_Count = max_targets_count,
    Min_Targets_miRNA = min_targets_mirna,
    Min_Targets_Count = min_targets_count
  ))
}

# 为所有基因组生成统计信息
generate_all_statistics <- function() {
  # 初始化结果数据框
  all_stats <- data.frame()
  
  # 为每个基因组计算统计信息
  for (genome in names(genome_data_list)) {
    if (!is.null(genome_data_list[[genome]])) {
      # 生成统计信息
      stats <- generate_genome_statistics(genome_data_list[[genome]], genome)
      
      # 添加到结果中
      all_stats <- rbind(all_stats, stats)
    }
  }
  
  # 保存统计信息到CSV文件
  write.csv(all_stats, "network_plots/network_statistics.csv", row.names = FALSE)
  
  # 打印统计信息摘要
  cat("\n基因组miRNA-靶基因网络统计信息:\n")
  print(all_stats)
  
  return(all_stats)
}
cat("生成网络统计信息...\n")
network_stats <- generate_all_statistics()
cat("网络统计信息已保存至 'network_plots/network_statistics.csv'\n")
cat("所有网络图保存完成。\n")