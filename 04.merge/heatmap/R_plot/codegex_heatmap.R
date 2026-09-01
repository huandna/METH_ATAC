# ============================
# 染色体信号热度图 - 优化颜色方案版
# 坐标：kb (1000bp) | 独立色阶 | 浅蓝到鲜红渐变 | 无空白
# 运行：Rscript heatmap.R 你的文件.csv
# ============================

library(ggplot2)
library(dplyr)
library(readr)
library(scales)  # 添加scales包，用于更好的颜色映射

# 接收外部文件参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("用法: Rscript heatmap.R 你的文件.csv")
}
file_path <- args[1]

# 读取数据
cat("📖 正在读取数据...\n")
data <- read_csv(file_path, show_col_types = FALSE)

# 数据清洗：保留Chr+数字染色体 + 坐标转 kb (除以1000)
cat("🔧 正在清洗数据...\n")
data_clean <- data %>%
  filter(grepl("^Chr\\d+$", chr)) %>%
  mutate(start_kb = start / 1000)  # 核心：bp → kb 转换

# 确保染色体按数字顺序排列
data_clean$chr_num <- as.numeric(gsub("Chr", "", data_clean$chr))
data_clean <- data_clean %>% arrange(chr_num, start_kb)
data_clean$chr <- factor(data_clean$chr, levels = unique(data_clean$chr))

# 检查并打印数据范围，用于调试
cat("📊 数据范围检查:\n")
cat(paste0("ATAC: ", min(data_clean$atac, na.rm=TRUE), " - ", max(data_clean$atac, na.rm=TRUE), "\n"))
cat(paste0("Methylation: ", min(data_clean$methyl, na.rm=TRUE), " - ", max(data_clean$methyl, na.rm=TRUE), "\n"))
cat(paste0("Activity Score: ", min(data_clean$activity_score, na.rm=TRUE), " - ", max(data_clean$activity_score, na.rm=TRUE), "\n"))

# 定义一个函数来绘制和保存热图，完成后立即释放内存
plot_and_save_heatmap <- function(data, signal_col, title, out_prefix) {
  cat(paste0("🎨 正在绘制", title, "热图...\n"))
  
# 获取数据范围
signal_values <- data[[signal_col]]
vmin <- min(signal_values, na.rm=TRUE)
vmax <- max(signal_values, na.rm=TRUE)

# 添加中位数和分位数检查
cat(paste0("  数据范围: ", vmin, " - ", vmax, "\n"))
cat(paste0("  中位数: ", median(signal_values, na.rm=TRUE), "\n"))
cat(paste0("  25%分位数: ", quantile(signal_values, 0.25, na.rm=TRUE), "\n"))
cat(paste0("  75%分位数: ", quantile(signal_values, 0.75, na.rm=TRUE), "\n"))
cat(paste0("  数据点数量: ", length(signal_values), "\n"))

# 如果数据点过多，进行抽样
if (nrow(data) > 50000) {
  cat("  ⚠️ 数据点过多，进行抽样...\n")
  data <- data %>% sample_frac(0.1)  # 抽取10%的数据
}

# 创建图形对象
p <- ggplot(data, aes(x = start_kb, y = chr)) +
  geom_tile(aes_string(fill = signal_col), height = 0.85, color = NA) +
  # 使用不含白色的颜色方案：深蓝 -> 中蓝 -> 粉红 -> 鲜红
  scale_fill_gradientn(colors = c("#004080", "#4080BF", "#FF6666", "#FF0000"), 
                       na.value = "gray90",
                       guide = guide_colorbar(title = title),
                       # 使用对数转换来处理数据范围
                       trans = "log1p") +
  labs(x = "Position (kb)", y = "Chromosome") +
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

  
  # 保存图片
  cat(paste0("💾 正在保存", title, "热图...\n"))
  ggsave(paste0(out_prefix, ".png"), p, width = 14, height = 5, dpi = 300)
  ggsave(paste0(out_prefix, ".pdf"), p, width = 14, height = 5)
  
  # 显式释放内存
  rm(p)
  gc()  # 强制垃圾回收
  
  cat(paste0("✅ 已保存: ", out_prefix, "\n"))
}

# 依次绘制和保存每个热图
plot_and_save_heatmap(data_clean, "atac", "ATAC", "ATAC_heatmap")
plot_and_save_heatmap(data_clean, "methyl", "Methylation", "Methylation_heatmap")
plot_and_save_heatmap(data_clean, "activity_score", "Activity Score", "Activity_score_heatmap")

cat("✅ 所有热图绘制完成！坐标：kb | 独立色阶 | 浅蓝到鲜红渐变\n")
