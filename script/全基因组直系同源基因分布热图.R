library(readxl)
library(tibble)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(dplyr)
library(tidyr)
library(cowplot)

## 设置工作目录和基因组顺序
setwd("~/Desktop/project/sRNA_prediction/比较基因组学解析鲜食和砧木葡萄群体间 miRNA 靶基因/")
order <- c("110R-hap1", "110R-hap2", "5BB-hap1", "5BB-hap2", "T2T", "MH-hap1", "MH-hap2", "SM-hap1", "SM-hap2")

## 读入数据
data <- read.csv("./data for visualization/Orthogroups_matrix.csv", sep = ";", header = TRUE, check.names = FALSE)
# data <- read.csv("../test/Orthogroups_matrix.csv", sep = ";", header = TRUE, check.names = FALSE)
message("原始数据维度: ", paste(dim(data), collapse = " x "))

## 剔除所有列（除第一列Orthogroup外）都为0的行
data <- data[rowSums(data[, -1] > 0) > 0, ]
message("剔除全为0的行后数据维度: ", paste(dim(data), collapse = " x "))

## 直系同源基因分布热图准备
color_function <- c("0" = "#a6cee3", "1" = "#fb9a99")
# 首先创建data_1变量，然后再设置行名
data_1 <- data[, -1]
rownames(data_1) <- data[, 1]
data_1 <- data_1[, order]  # 按指定顺序排列列
data_1[] <- lapply(data_1, function(x) as.numeric(as.character(x)))
data_1 <- as.matrix(data_1)
data_1[data_1 > 0] <- 1  # 将所有大于0的值设为1

## 绘制热图
p1 <- Heatmap(
  data_1,
  heatmap_legend_param = list(
    at = c(1, 0),
    labels = c("Present", "Absent"),
    title = "",
    legend_direction = "horizontal",
    nrow = 1,
    ncol = 2
  ),
  col = color_function,
  show_row_names = FALSE,
  border = FALSE,
  row_names_gp = gpar(fontsize = 3),
  column_names_gp = gpar(fontsize = 10, just = "centre"),
  column_names_centered = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "bottom",
  column_names_rot = 0
)
## 单独保存p1为PNG格式（图例在顶部，水平排列）
png("Fig.2_heatmap.png", width = 8.26, height = 4, units = "in", res = 600)
draw(p1, 
     heatmap_legend_side = "top",    # 图例位置在顶部
     annotation_legend_side = "top"  # 注释图例也在顶部
)
dev.off()

## 计算核心基因数量（在所有物种中都至少有一个基因的 OG）
core_genes <- data %>%
  rowwise() %>%
  filter(all(c_across(-Orthogroup) > 0)) %>%  # 所有物种基因数都大于0
  nrow()
message("核心基因数量: ", core_genes)

## 计算各种基因统计数据

# 每个物种总基因数量
total_genes_per_species <- colSums(data[-1]) 

# 每个物种的 PAV 基因数量
pav_genes_per_species <- colSums(data[-1] == 0)  

## 计算每个物种的特有基因数量
unique_genes_per_species <- sapply(
  colnames(data)[-1],
  function(species) {
    sum(apply(data[-1], 1, function(row) {
      row[species == colnames(data)[-1]] > 0 & 
        sum(row[-which(colnames(data)[-1] == species)]) == 0
    }))
  }
)

## 计算百分比
pav_percentage <- pav_genes_per_species / total_genes_per_species * 100
unique_percentage <- unique_genes_per_species / total_genes_per_species * 100

## 整理统计结果
result <- data.frame(
  Species = colnames(data)[-1],
  Core_Genes = core_genes,
  PAV_Genes = pav_genes_per_species,
  Unique_Genes = unique_genes_per_species,
  Total_Genes = total_genes_per_species,
  PAV_Percentage = pav_percentage,
  Unique_Percentage = unique_percentage
)

## 转换为长格式用于绘图
plot_data <- result %>%
  select(Species, Core_Genes, PAV_Genes, Unique_Genes) %>%
  pivot_longer(
    cols = -Species,
    names_to = "Category",
    values_to = "Count"
  ) %>%
  group_by(Species) %>%
  mutate(Percentage = Count / sum(Count))

## 指定三种基因的顺序
plot_data <- plot_data %>%
  mutate(Category = factor(Category, levels = rev(c("Core_Genes", "PAV_Genes", "Unique_Genes"))))

## 指定基因组顺序
plot_data$Species <- factor(plot_data$Species, levels = order)

## 绘制柱状图
p2 <- ggplot(plot_data, aes(x = Species, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3) +
  labs(
    title = "",
    x = "",
    y = "",
    fill = ""
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(
    values = rev(c("#ccebc5", "#fed9a6", "#b3cde3")),
    labels = rev(c("Core genes", "PAV genes", "Specific genes"))
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
    axis.ticks.x = element_line(color = "black", size = .5),
    axis.ticks.y = element_line(color = "black", size = .5),
    axis.line.x = element_line(color = "black", size = .5)
  )
## 单独保存p2为PNG格式
png("Fig.2_barplot.png", width = 8.26, height = 3, units = "in", res = 600)
print(p2)
dev.off()

## 合并图像
p1_grob <- grid.grabExpr(draw(p1, heatmap_legend_side = "top", annotation_legend_side = "top"))
p <- plot_grid(
  ggdraw(p1_grob),
  p2,
  ncol = 1,
  nrow = 2,
  align = "v",
  rel_heights = c(2, 1)
)
png("PAV.png", width = 8.26, height = 6, units = "in", res = 600)
print(p1_grob)
dev.off()
## 保存图像
# png("Fig.2.png", width = 9.5, height = 8, units = "in", res = 600)
# print(p)
# dev.off()
