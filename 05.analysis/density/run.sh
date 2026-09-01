# 统计第4列（ATAC信号）的分布
# python plot_signal_distribution.py --input genome_ml_dataset.bed --column 4 --signal_type ATAC

# 统计第5列（甲基化信号）的分布
#python plot_signal_distribution.py --input genome_ml_dataset.bed --column 5 --signal_type Methylation
#
python plot_signal_distribution.py --input ./atac_matrix.bed  --signal_type atac --matrix_format

python plot_signal_distribution.py --input ./methyl_matrix.bed  --signal_type meth --matrix_format

