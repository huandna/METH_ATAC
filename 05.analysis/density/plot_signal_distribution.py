import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rcParams

# 设置中文字体
rcParams['font.sans-serif'] = ['DejaVu Sans']
rcParams['axes.unicode_minus'] = False

def load_data(file_path, signal_col, selected_chroms=None, min_signal=None, is_matrix_format=False):
    """加载数据文件，支持两种格式：
    1. 标准格式：每行一个位置，第一列为染色体ID，指定列为信号值
    2. 矩阵格式：每行一个染色体，第一列为染色体ID，其余列为信号值
    
    参数:
        file_path: 文件路径
        signal_col: 要统计的信号列（从1开始计数，仅用于标准格式）
        selected_chroms: 仅显示指定的序列ID列表
        min_signal: 仅统计信号值大于等于此值的记录
        is_matrix_format: 是否为矩阵格式
    """
    data = []
    
    if is_matrix_format:
        # 矩阵格式：每行一个染色体，第一列为染色体ID，其余列为信号值
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                
                chrom = parts[0]
                
                # 检查染色体ID格式：chr或Chr开头，后跟数字
                if not (chrom.startswith('chr') or chrom.startswith('Chr')):
                    continue
                
                # 提取数字部分并检查是否为有效数字
                chrom_num = chrom[3:] if chrom.startswith('Chr') else chrom[3:]
                if not chrom_num.isdigit():
                    continue
                
                # 如果指定了要选择的序列ID，则只保留这些序列ID
                if selected_chroms is not None and chrom not in selected_chroms:
                    continue
                
                # 提取所有信号值
                signals = [float(val) for val in parts[1:]]
                
                # 如果指定了最小信号值，则过滤掉小于该值的记录
                if min_signal is not None:
                    signals = [s for s in signals if s >= min_signal]
                    # 如果过滤后没有信号值，跳过该染色体
                    if not signals:
                        continue
                
                # 将每个信号值作为一个记录
                for signal in signals:
                    data.append([chrom, signal])
    else:
        # 标准格式：每行一个位置，第一列为染色体ID，指定列为信号值
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= signal_col:
                    chrom = parts[0]
                    
                    # 如果指定了要选择的序列ID，则只保留这些序列ID
                    if selected_chroms is not None:
                        if chrom not in selected_chroms:
                            continue
                    else:
                        # 否则，检查染色体ID格式：chr或Chr开头，后跟数字
                        if not (chrom.startswith('chr') or chrom.startswith('Chr')):
                            continue
                        
                        # 提取数字部分并检查是否为有效数字
                        chrom_num = chrom[3:] if chrom.startswith('Chr') else chrom[3:]
                        if not chrom_num.isdigit():
                            continue
                    
                    signal = float(parts[signal_col-1])  # Python索引从0开始，所以需要减1
                    
                    # 如果指定了最小信号值，则过滤掉小于该值的记录
                    if min_signal is not None and signal < min_signal:
                        continue
                    
                    data.append([chrom, signal])
    
    df = pd.DataFrame(data, columns=['chrom', 'signal'])
    return df




def plot_signal_distribution(df, signal_type, output_dir):
    """绘制信号分布曲线"""
    # 获取所有染色体
    chromosomes = df['chrom'].unique()
    chromosomes.sort()
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 为每条染色体绘制分布曲线
    colors = plt.cm.tab20(np.linspace(0, 1, len(chromosomes)))
    
    for i, chrom in enumerate(chromosomes):
        chrom_data = df[df['chrom'] == chrom]['signal']
        # 计算核密度估计
        if len(chrom_data) > 0 and chrom_data.sum() > 0:  # 只处理有数据的染色体
            try:
                density = np.histogram(chrom_data, bins=100, density=True)[0]
                bins = np.histogram(chrom_data, bins=100, density=True)[1]
                bin_centers = (bins[:-1] + bins[1:]) / 2
                
                ax.plot(bin_centers, density, label=chrom, color=colors[i], linewidth=2)
            except Exception as e:
                print(f"处理染色体 {chrom} 时出错: {e}")
    
    # 设置图形属性
    ax.set_xlabel(f'{signal_type} Signal Value', fontsize=14, fontweight='bold')
    ax.set_ylabel('Density', fontsize=14, fontweight='bold')
    ax.set_title(f'{signal_type} Signal Distribution by Chromosome', fontsize=16, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.7)
    
    # 保存图形
    os.makedirs(output_dir, exist_ok=True)
    pdf_path = os.path.join(output_dir, f'{signal_type}_distribution.pdf')
    png_path = os.path.join(output_dir, f'{signal_type}_distribution.png')
    
    plt.savefig(pdf_path, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(png_path, format='png', dpi=300, bbox_inches='tight')
    
    print(f"✅ 已保存: {pdf_path}")
    print(f"✅ 已保存: {png_path}")
    
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='统计染色体信号值频率分布并绘制曲线')
    parser.add_argument('--input', required=True, help='输入文件路径')
    parser.add_argument('--column', type=int, help='要统计的信号列（从1开始计数，仅用于标准格式）')
    parser.add_argument('--signal_type', required=True, help='信号类型名称（用于标题和文件名）')
    parser.add_argument('--out_dir', default='signal_distribution', help='输出目录')
    parser.add_argument('--selected_chroms', nargs='+', help='仅显示指定的序列ID（如 Chr01 Chr02 Chr03）')
    parser.add_argument('--min_signal', type=float, help='仅统计信号值大于等于此值的记录')
    parser.add_argument('--matrix_format', action='store_true', help='输入文件为矩阵格式（每行一个染色体，第一列为染色体ID，其余列为信号值）')
    args = parser.parse_args()
    
    # 检查参数
    if not args.matrix_format and args.column is None:
        parser.error("--column 参数在非矩阵格式下是必需的")
    
    # 加载数据
    print(f"📖 正在读取数据: {args.input}")
    print(f"📝 文件格式: {'矩阵格式' if args.matrix_format else '标准格式'}")
    df = load_data(args.input, args.column, args.selected_chroms, args.min_signal, args.matrix_format)
    
    # 打印数据统计信息
    print(f"📊 数据统计:")
    print(f"  总记录数: {len(df)}")
    print(f"  染色体数量: {df['chrom'].nunique()}")
    if len(df) > 0:
        print(f"  信号范围: {df['signal'].min():.4f} - {df['signal'].max():.4f}")
        print(f"  信号均值: {df['signal'].mean():.4f}")
        print(f"  信号中位数: {df['signal'].median():.4f}")
    else:
        print("  ⚠️ 警告: 没有符合条件的数据！")
        return
    
    # 如果指定了选择的序列ID，打印出来
    if args.selected_chroms is not None:
        print(f"  选择的序列ID: {', '.join(args.selected_chroms)}")
    
    # 如果指定了最小信号值，打印出来
    if args.min_signal is not None:
        print(f"  最小信号值阈值: {args.min_signal}")
    
    # 绘制分布曲线
    print(f"🎨 正在绘制分布曲线...")
    plot_signal_distribution(df, args.signal_type, args.out_dir)
    
    print(f"✅ 完成！分布曲线已保存到 {args.out_dir}")

if __name__ == '__main__':
    main()


