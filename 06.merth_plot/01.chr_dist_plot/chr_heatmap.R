#!/usr/bin/env Rscript

# 加载必要的包
library(ggplot2)
library(reshape2)
library(dendextend)
library(ggdendro)
library(grid)

# 读取数据
data <- read.csv("methylation_level_by_chromosome.csv", row.names = 1)

# 处理重复样本：取平均值
# 假设样本名格式为 "sample_rep" 或类似，这里简单取唯一样本名
# 如果样本名已经唯一，可以跳过这一步
unique_samples <- unique(rownames(data))
if (length(unique_samples) < nrow(data)) {
  # 如果有重复样本，计算平均值
  averaged_data <- aggregate(. ~ rownames(data), data = data, FUN = mean)
  rownames(averaged_data) <- averaged_data[,1]
  averaged_data <- averaged_data[,-1]
  data <- averaged_data
}

# 计算样本相关性
cor_matrix <- cor(t(data), method = "pearson")

# 创建层次聚类树
hc <- hclust(dist(1 - cor_matrix), method = "complete")

# 创建树状图数据
dend <- as.dendrogram(hc)
dend_data <- dendro_data(dend)

# 按聚类顺序重新排列数据
ordered_data <- data[labels(hc), ]

# 将数据转换为长格式，用于ggplot2
melted_data <- melt(as.matrix(ordered_data))
colnames(melted_data) <- c("Sample", "Chromosome", "Methylation_Level")

# 创建热图
p1 <- ggplot(melted_data, aes(x = Chromosome, y = Sample, fill = Methylation_Level)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0.5, limits = c(0, 1),
                       name = "Methylation\nLevel") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        panel.grid = element_blank())

# 创建树状图
p2 <- ggplot(segment(dend_data)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_reverse(expand = c(0.2, 0)) +
  theme_dendro()

# 合并树状图和热图
library(gridExtra)
final_plot <- grid.arrange(p2, p1, ncol = 2, widths = c(1, 4))

# 保存图片
ggsave("methylation_heatmap_with_dendrogram.pdf", final_plot, width = 12, height = 10)

print("热图已保存到 methylation_heatmap_with_dendrogram.pdf")
