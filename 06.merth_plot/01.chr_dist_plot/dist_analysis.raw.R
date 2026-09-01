#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)

# 读取数据
data <- read.csv("methylation_distribution_all.csv")

# 提取样本基础名称（去除末尾的数字）
melted_data <- melt(data, id.vars = c("Methylation_Level", "Context_Type"), 
                    variable.name = "Sample", value.name = "Frequency")

# 从样本名中提取基础名称（去除末尾的数字）
melted_data$Base_Sample <- gsub("\\d+$", "", melted_data$Sample)

# 按基础样本名、上下文类型和甲基化水平分组，计算频率的平均值
averaged_data <- melted_data %>%
  group_by(Base_Sample, Context_Type, Methylation_Level) %>%
  summarise(Frequency = mean(Frequency, na.rm = TRUE), .groups = "drop")

# 确保甲基化水平按顺序排列
averaged_data$Methylation_Level <- factor(averaged_data$Methylation_Level, 
                                          levels = c("0", "10", "20", "30", "40", "50", 
                                                    "60", "70", "80", "90", "100"))

# 获取样本数量
num_samples <- length(unique(averaged_data$Base_Sample))

# 选择配色方案 - 使用Set2或Set3配色，这些配色对比度较高且颜色区分明显
if (num_samples <= 8) {
  colors <- brewer.pal(num_samples, "Set2")
} else if (num_samples <= 12) {
  colors <- brewer.pal(num_samples, "Set3")
} else {
  # 如果样本数超过12，使用Paired配色并循环使用
  colors <- rep(brewer.pal(12, "Paired"), length.out = num_samples)
}

# 创建折线图，按上下文类型分面
p <- ggplot(averaged_data, aes(x = Methylation_Level, y = Frequency, 
                               group = Base_Sample, color = Base_Sample)) +
  geom_line(size = 0.8, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~ Context_Type, nrow = 1, scales = "free_y") +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(values = colors) +
  labs(title = "Methylation Level Distribution by Context Type",
       x = "Methylation Level (%)",
       y = "Frequency (%)",
       color = "Sample") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  )

# 保存为PDF
ggsave("methylation_distribution.pdf", plot = p, width = 18, height = 6, device = "pdf")

print("折线图已保存到 methylation_distribution.pdf")
