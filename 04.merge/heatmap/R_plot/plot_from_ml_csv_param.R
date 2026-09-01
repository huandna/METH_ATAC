# ============================
# 染色体信号热度图 - 最终定稿版
# 坐标：kb (1000bp) | 独立色阶 | 蓝白红渐变 | 无空白
# 运行：Rscript heatmap.R 你的文件.csv
# ============================

library(ggplot2)
library(dplyr)
library(readr)

# 接收外部文件参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("用法: Rscript heatmap.R 你的文件.csv")
}
file_path <- args[1]

# 读取数据
data <- read_csv(file_path, show_col_types = FALSE)

# 数据清洗：保留Chr+数字染色体 + 坐标转 kb (除以1000)
data_clean <- data %>%
  filter(grepl("^Chr\\d+$", chr)) %>%
  mutate(start_kb = start / 1000)  # 核心：bp → kb 转换
data_clean$chr <- factor(data_clean$chr)

# ================== 1. ATAC 独立色阶绘图 ==================
p1 <- ggplot(data_clean, aes(x = start_kb, y = chr)) +
  geom_tile(aes(fill = atac), height = 0.85, color = NA) +
  scale_fill_gradientn(colors = c("#0055FF", "white", "#FF2222"), na.value = "gray90") +
  labs(x = "Position (kb)", y = "Chromosome", fill = "ATAC") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(linewidth = 0.8),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# ================== 2. 甲基化 独立色阶绘图 ==================
p2 <- ggplot(data_clean, aes(x = start_kb, y = chr)) +
  geom_tile(aes(fill = methyl), height = 0.85, color = NA) +
  scale_fill_gradientn(colors = c("#0055FF", "white", "#FF2222"), na.value = "gray90") +
  labs(x = "Position (kb)", y = "Chromosome", fill = "Methylation") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(linewidth = 0.8),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# ================== 3. 活性分数 独立色阶绘图 ==================
p3 <- ggplot(data_clean, aes(x = start_kb, y = chr)) +
  geom_tile(aes(fill = activity_score), height = 0.85, color = NA) +
  scale_fill_gradientn(colors = c("#0055FF", "white", "#FF2222"), na.value = "gray90") +
  labs(x = "Position (kb)", y = "Chromosome", fill = "Activity Score") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(linewidth = 0.8),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# 保存图片：PNG + PDF
ggsave("ATAC_heatmap.png", p1, width = 14, height = 5, dpi = 300)
ggsave("ATAC_heatmap.pdf", p1, width = 14, height = 5)
ggsave("Methylation_heatmap.png", p2, width = 14, height = 5, dpi = 300)
ggsave("Methylation_heatmap.pdf", p2, width = 14, height = 5)
ggsave("Activity_score_heatmap.png", p3, width = 14, height = 5, dpi = 300)
ggsave("Activity_score_heatmap.pdf", p3, width = 14, height = 5)

cat("✅ 绘图完成！坐标：kb | 独立色阶 | 蓝白红渐变\n")