#!/usr/bin/env Rscript

# 加载必要的包
library(ggplot2)
library(reshape2)

# 读取数据
data <- read.csv("methylation_count_by_chromosome_averaged.csv", row.names = 1)

# 转置数据，使染色体为行，样本为列
data_t <- as.data.frame(t(data))

# 添加染色体列
data_t$Chromosome <- rownames(data_t)

# 提取染色体编号并排序
data_t$ChrNum <- as.numeric(gsub("Chr", "", data_t$Chromosome))
data_t <- data_t[order(data_t$ChrNum), ]

# 将数据转换为长格式，用于ggplot2
melted_data <- melt(data_t, id.vars = c("Chromosome", "ChrNum"), 
                    variable.name = "Sample", value.name = "Count")

# 创建折线图
# 创建折线图
# p <- ggplot(melted_data, aes(x = Chromosome, y = Count, group = Sample, color = Sample)) +
#   geom_line(size = 1) +
#   geom_point(size = 2) +
#   labs(title = "Number of Methylation Sites (>0.5) by Chromosome",
#        x = "Chromosome",
#        y = "Count of Methylation Sites",
#        color = "Sample") +
#   theme_minimal() +
#   theme(
#     # 移除背景网格线
#     panel.grid = element_blank(),
#     # 添加左边和下边的坐标轴
#     axis.line = element_line(color = "black", linewidth = 0.5),
#     # 移除顶部和右侧边框
#     axis.line.x.top = element_blank(),
#     axis.line.y.right = element_blank(),
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "right",
#     plot.title = element_text(hjust = 0.5, size = 14),
#     axis.title = element_text(size = 12),
#     legend.title = element_text(size = 12)
#   )

# # 创建柱状图（dodge）
# p <- ggplot(melted_data, aes(x = Chromosome, y = Count, fill = Sample)) +
#   geom_col(position = position_dodge(width = 0.7), width = 0.5) +
#   labs(title = "Number of Methylation Sites (>0.5) by Chromosome",
#        x = "Chromosome",
#        y = "Count of Methylation Sites",
#        fill = "Sample") +
#   theme_minimal() +
#   scale_fill_brewer(palette = "Set1") +  # 使用Set1配色方案，颜色区分度更强
#   theme(
#     # 移除背景网格线
#     panel.grid = element_blank(),
#     # 添加左边和下边的坐标轴
#     axis.line = element_line(color = "black", linewidth = 0.5),
#     # 移除顶部和右侧边框
#     axis.line.x.top = element_blank(),
#     axis.line.y.right = element_blank(),
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "right",
#     plot.title = element_text(hjust = 0.5, size = 14),
#     axis.title = element_text(size = 12),
#     legend.title = element_text(size = 12)
#   )
# 创建柱状图（stack）
# 创建柱状图
#scale_fill_manual(values = c("#3C5488", "#00A087", "#F39B7F", "#E64B35"))



# 创建柱状图
# 创建柱状图
p <- ggplot(melted_data, aes(x = Chromosome, y = Count, fill = Sample)) +
  geom_col(position = "stack", width = 0.5, color = "black", size = 0.2) +  # 添加黑色边框，边框宽度为0.2
  labs(title = "Number of Methylation Sites (>0.5) by Chromosome",
       x = "Chromosome",
       y = "Count of Methylation Sites (×10⁴)",  # 添加纵坐标单位
       fill = "Sample") +
  theme_minimal(base_size = 8) +  # 设置基础字体大小为8，使整体变小
  scale_fill_manual(values = c("#3C5488", "#00A087", "#F39B7F", "#E64B35")) +
  scale_y_continuous(labels = function(x) x/10000) +  # 将纵坐标数值除以10000
  theme(
    # 移除背景网格线
    panel.grid = element_blank(),
    # 添加左边和下边的坐标轴
    axis.line = element_line(color = "black", linewidth = 0.3),
    # 移除顶部和右侧边框
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),  # x轴文字变大变粗
    axis.text.y = element_text(size = 10, face = "bold"),  # y轴文字变大变粗
    axis.title.x = element_text(size = 11, face = "bold"),  # x轴标题变大变粗
    axis.title.y = element_text(size = 11, face = "bold"),  # y轴标题变大变粗
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )



# 保存为PDF
ggsave("methylation_stats_bar.pdf", plot = p, width = 12, height = 6, device = "pdf")

print("折线图已保存到 methylation_stats_bar.pdf")
