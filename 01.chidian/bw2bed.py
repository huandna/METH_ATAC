import os
import glob
import pyBigWig
import numpy as np
import argparse

def bw_to_bed_with_signal(bw_path, output_dir, bin_size=10):
    """
    BW转BED，强制保留信号值（第四列）
    :param bw_path: 输入BW文件路径
    :param output_dir: 输出目录
    :param bin_size: 分辨率（bp）
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 提取文件名，生成输出BED路径
    bw_basename = os.path.basename(bw_path)
    bed_basename = bw_basename.replace(".bw", ".bed")
    bed_path = os.path.join(output_dir, bed_basename)
    
    # 打开BW文件
    try:
        bw = pyBigWig.open(bw_path)
        if bw is None:
            print(f"❌ 无法打开文件：{bw_path}")
            return
        print(f"✅ 开始转换：{bw_basename} → {bed_basename}")
    except Exception as e:
        print(f"❌ 打开文件失败：{e}")
        return
    
    # 获取BW文件的染色体和长度（优先用BW内置的，避免依赖外部chrom.sizes）
    chrom_sizes = bw.chroms()
    if not chrom_sizes:
        print(f"❌ {bw_basename} 无染色体信息，跳过")
        bw.close()
        return
    
    # 写入BED文件（强制保留第四列信号值）
    with open(bed_path, 'w') as bed_f:
        for chrom, chrom_size in chrom_sizes.items():
            # 仅跳过明显的非标准染色体（避免无意义数据，可注释掉这行）
            if 'random' in chrom or 'chrM' in chrom or 'Un' in chrom:
                continue
            
            # 按bin_size遍历每个区间，无任何过滤
            for start in range(0, chrom_size, bin_size):
                end = min(start + bin_size, chrom_size)
                if start >= end:
                    continue
                
                # 提取区间信号，无信号则赋值为0.0
                sig_list = bw.values(chrom, start, end)
                if sig_list is None or len(sig_list) == 0:
                    avg_signal = 0.0
                else:
                    # 计算均值，NaN/Inf全部替换为0.0
                    sig_array = np.array(sig_list, dtype=np.float64)
                    sig_array = np.nan_to_num(sig_array, nan=0.0, posinf=0.0, neginf=0.0)
                    avg_signal = np.mean(sig_array)
                
                # 强制输出4列BED：chrom start end signal（信号值必存在）
                bed_f.write(f"{chrom}\t{start}\t{end}\t{avg_signal:.4f}\n")
    
    # 关闭文件
    bw.close()
    print(f"✅ 转换完成：{bed_path}（每一行都有信号值）")

def main():
    # 简单的命令行参数（新手友好）
    parser = argparse.ArgumentParser(description='BW转BED（强制保留信号值）')
    parser.add_argument('--input', '-i', required=True, help='BW文件路径/目录（如 ./test.bw 或 ./bw_dir）')
    parser.add_argument('--output', '-o', default='./bw2bed_results', help='输出目录')
    parser.add_argument('--bin_size', '-b', type=int, default=10, help='分辨率（bp）')
    
    args = parser.parse_args()
    
    # 处理输入：单个文件 or 目录下所有BW文件
    if os.path.isfile(args.input) and args.input.endswith('.bw'):
        # 转换单个BW文件
        bw_to_bed_with_signal(args.input, args.output, args.bin_size)
    elif os.path.isdir(args.input):
        # 批量转换目录下所有BW文件
        bw_files = glob.glob(os.path.join(args.input, "*.bw"))
        if len(bw_files) == 0:
            print(f"❌ 目录 {args.input} 下未找到BW文件")
            return
        print(f"📁 找到 {len(bw_files)} 个BW文件，开始批量转换...")
        for bw_file in bw_files:
            bw_to_bed_with_signal(bw_file, args.output, args.bin_size)
    else:
        print(f"❌ 输入无效：{args.input}（需为BW文件或目录）")

if __name__ == "__main__":
    main()