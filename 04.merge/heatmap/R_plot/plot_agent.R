# 染色体信号值热图绘制脚本
# 使用data.table和ggplot2高效处理大规模数据

# 加载必要的包
library(data.table)
library(ggplot2)

# 定义颜色渐变函数 - 从浅蓝到鲜红
get_color_gradient <- function() {
  # 使用自定义颜色渐变：浅蓝(#ADD8E6) -> 中蓝(#4169E1) -> 橙色(#FF8C00) -> 鲜红(#FF0000)
  colors <- c("#ADD8E6", "#4169E1", "#FF8C00", "#FF0000")
  return(colors)
}

# 处理数据并绘制热图的函数
process_and_plot <- function(input_file, output_prefix) {
  # 使用data.table高效读取数据
  message("正在读取数据...")
  dt <- fread(input_file, header = TRUE, sep = ",", showProgress = TRUE)
  
  # 获取列名
  col_names <- colnames(dt)
  
  # 假设前两列是染色体和位置，最后三列是信号值
  chr_col <- col_names[1]  # 染色体列
  pos_col <- col_names[2]  # 位置列
  signal_cols <- col_names[(length(col_names)-2):length(col_names)]  # 最后三列信号值
  
  # 将位置转换为kb单位
  message("正在转换位置单位...")
  dt[, (pos_col) := get(pos_col) / 1000]
  
  # 确保染色体是因子型并按顺序排列
  message("正在处理染色体信息...")
  dt[, (chr_col) := as.factor(get(chr_col))]
  
  # 获取所有染色体并按顺序排列
  chromosomes <- sort(unique(dt[[chr_col]]))
  dt[, (chr_col) := factor(get(chr_col), levels = chromosomes)]
  
  # 为每个信号列创建热图
  for (i in 1:length(signal_cols)) {
    signal_col <- signal_cols[i]
    message(paste("正在处理信号列:", signal_col))
    
    # 获取当前信号列的范围
    min_val <- min(dt[[signal_col]], na.rm = TRUE)
    max_val <- max(dt[[signal_col]], na.rm = TRUE)
    message(paste("信号值范围:", min_val, "到", max_val))
    
    # 创建ggplot热图
    p <- ggplot(dt, aes_string(x = pos_col, y = chr_col, fill = signal_col)) +
      geom_tile() +
      scale_fill_gradientn(colors = get_color_gradient(), 
                          limits = c(min_val, max_val),
                          na.value = "grey50") +
      labs(x = "位置 (kb)", y = "染色体", fill = signal_col) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.grid.major = element_blank(),
        legend.grid.minor = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)
      ) +
      ggtitle(paste("染色体信号热图 -", signal_col))
    
    # 保存热图
    output_file <- paste0(output_prefix, "_", signal_col, ".png")
    message(paste("正在保存热图到:", output_file))
    ggsave(filename = output_file, plot = p, width = 12, height = 10, dpi = 300)
    
    # 释放内存
    rm(p)
    gc()
  }
  
  # 释放数据内存
  rm(dt)
  gc()
  message("处理完成!")
}

# 从命令行参数获取输入文件和输出前缀
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("使用方法: Rscript plot_chromosome_heatmap.R <input_file> [output_prefix]")
}

input_file <- args[1]
output_prefix <- if (length(args) >= 2) args[2] else "chromosome_heatmap"

# 执行主函数
process_and_plot(input_file, output_prefix)
