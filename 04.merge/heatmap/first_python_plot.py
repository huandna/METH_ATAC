import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ------------------------------------------------------------------------------
# 读取矩阵
# ------------------------------------------------------------------------------
def load_full_matrix(file_path):
    chroms = []
    data_list = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            chrom = parts[0]
            vals = np.array(parts[1:], dtype=np.float32)
            chroms.append(chrom)
            data_list.append(vals)
    return chroms, data_list
def plot_strip_heatmap(chroms, data_list, title, out_path, bin_size=10, cmap_type='default'):
    max_len = max(len(d) for d in data_list)
    num_chroms = len(chroms)
    all_vals = np.concatenate(data_list)

    # 自动计算颜色范围（2%~98% 极值截断）
    v_low = np.percentile(all_vals, 2)
    v_high = np.percentile(all_vals, 98)
    v_mid = (v_low + v_high) / 2.0

    # 自定义配色：蓝 → 白 → 红（中间白，区分正负最清晰）
    if cmap_type == 'divergent':
        colors = ['blue', 'white', 'red']
        nodes = [0.0, 0.5, 1.0]
        cmap = LinearSegmentedColormap.from_list('custom_divergent', list(zip(nodes, colors)))
        vmin, vmax = v_low, v_high
    else:
        cmap = 'Reds' if 'ATAC' in title else 'Blues'
        vmin, vmax = v_low, v_high

    # 画布大小
    fig, axes = plt.subplots(num_chroms, 1, figsize=(28, num_chroms * 1.0), sharex=True)
    if num_chroms == 1:
        axes = [axes]

    # 绘图
    for ax, chrom, data in zip(axes, chroms, data_list):
        ax.set_yticks([])
        ax.set_ylabel(chrom, fontsize=16, weight='bold', rotation=0, ha='right', va='center')
        im = ax.imshow(
            data[np.newaxis, :],
            aspect='auto',
            cmap=cmap,
            vmin=vmin,
            vmax=vmax
        )
        ax.set_xlim(0, max_len)
        ax.tick_params(axis='x', labelsize=14)
        # 去除边框
        for spine in ax.spines.values():
            spine.set_visible(False)

    # 底部X轴：转 kb，固定为整数坐标
    # 底部X轴：转 kb，固定为整数坐标，只显示几个关键刻度
# 底部X轴：转 mb，固定为整数坐标，只显示几个关键刻度
    max_mb = int(max_len * bin_size / 1000000)
    # 计算合适的刻度间隔，确保只显示5-6个刻度
    num_ticks = 6
    step = max(1, max_mb // (num_ticks - 1))
    # 确保step是1或0.5的倍数，便于阅读
    if step >= 1:
        step = round(step)
    else:
        step = 0.5

    x_ticks = np.arange(0, max_mb + step, step)
    # 转换回数据索引
    x_indices = x_ticks * 1000000 / bin_size
    x_labels = [f"{x} mb" for x in x_ticks]

    axes[-1].set_xticks(x_indices)
    axes[-1].set_xticklabels(x_labels, fontsize=15, weight='bold')
    # 只保留底部坐标轴
    axes[-1].spines['bottom'].set_visible(True)
    axes[-1].spines['left'].set_visible(False)
    axes[-1].spines['top'].set_visible(False)
    axes[-1].spines['right'].set_visible(False)


    # 总标题
    fig.suptitle(title, fontsize=22, weight='bold', y=0.98)

    # 右上角添加图例（颜色条），缩短长度
# 右上角添加图例（颜色条），垂直方向，从上到下渐变
    cbar_ax = fig.add_axes([0.85, 0.85, 0.03, 0.1])
    cb = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
    cb.ax.tick_params(labelsize=12)
    cb.set_label('Signal', fontsize=14, weight='bold')

    plt.tight_layout()
    plt.subplots_adjust(top=0.94)
    plt.savefig(out_path, dpi=300, format='pdf', bbox_inches='tight')
    plt.close()
    print(f"✅ 保存：{out_path}")



# ------------------------------------------------------------------------------
# 主程序
# ------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description='独立染色体水平条热图（最终PDF版）')
    parser.add_argument('--atac', required=True, help='atac_matrix.bed')
    parser.add_argument('--methyl', required=True, help='methyl_matrix.bed')
    parser.add_argument('--score', required=True, help='activity_score_matrix.bed')
    parser.add_argument('--out_dir', default='genome_pdf_heatmaps', help='输出目录')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    chroms, atac_data = load_full_matrix(args.atac)
    _, methyl_data = load_full_matrix(args.methyl)
    _, score_data = load_full_matrix(args.score)

    # 分别绘制三张独立图
    plot_strip_heatmap(chroms, atac_data,  'ATAC Signal',
        os.path.join(args.out_dir, 'atac_heatmap.pdf'), cmap_type='reds')
    plot_strip_heatmap(chroms, methyl_data, 'Methylation Level',
        os.path.join(args.out_dir, 'methyl_heatmap.pdf'), cmap_type='blues')
    plot_strip_heatmap(chroms, score_data,  'Activity Score (ATAC*(1-Methyl))',
        os.path.join(args.out_dir, 'activity_score_heatmap.pdf'), cmap_type='divergent')

    print("\n🎉 全部完成！3张PDF高清热图已生成！")

if __name__ == '__main__':
    main()