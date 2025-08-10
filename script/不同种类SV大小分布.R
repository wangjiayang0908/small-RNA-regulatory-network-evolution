library(ggplot2)
library(readxl)
library(tidyverse)
library(dplyr)
library(viridisLite)
library(viridis)
library(hrbrthemes)
library(ggbeeswarm)
library(ggpubr)
library(tidyverse)

setwd("~/Desktop/project/sRNA_prediction/比较基因组学解析鲜食和砧木葡萄群体间 miRNA 靶基因/")
box_data_2 <- read.csv("./data for visualization/merged_sv_distribution.csv",
                       header = T,sep = ',',check.names = F)
box_data_2 <- gather(box_data_2,len,count,-type)
# head(box_data_2)

p2 <-ggplot(box_data_2, aes(x=len,y=count,fill=type)) + 
  geom_bar(stat = "identity", position = "dodge") +
  # ylim(0,2000)+
  # scale_fill_manual(values=light_colors) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.7),  # 和柱形图对齐
    vjust = -0.2,  # 标签的位置，相对于柱子的顶部
    size = 3       # 标签字体大小
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    # panel.grid.major.x = element_line(color = "black",size = .2,linetype = "dashed"),
    axis.line.x=element_blank(),
    axis.line.y=element_blank(),
    axis.text.x=element_text(angle = 45,hjust = 1,vjust=1,size = 7,color = "black"),
    axis.text.y = element_text(hjust = 1,vjust=.5,size = 7,color = "black"),
    axis.ticks.x = element_line(color = "black",size = .1),
    axis.ticks.y = element_line(color = "black",size = .1),
    strip.text = element_text(size = 8,color = "black",face = "bold",hjust = 0.5),#调整每个分组标题的字体大小
    strip.background = element_blank(), # 隐藏子图标题部分的灰色背景
    panel.border = element_rect(color = "black", fill = NA, size = .3), # 增加边框
  #   legend.position = "top"
  # )+
  legend.position = "none") +
  guides(fill=guide_legend(),color="none")+
  # coord_flip() +
  labs(
    x = "",
    y = "",
    title = ""
  ) +
  facet_wrap(~type, nrow  = 3,scales = "free_y")+
  scale_x_discrete(limits = c("50-100","100-500", "500-1k", "1k-5k",
                                  "5k-10k", "10k-100k", "100k-1M", ">1M"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)));p2
png("Fig.1.png", width = 3, height = 9, units = "in", res = 600)
print(p2) 
dev.off()
