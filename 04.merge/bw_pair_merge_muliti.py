import os
import sys
import time
import glob
import argparse
import numpy as np
import pyBigWig
import subprocess
from collections import defaultdict
from datetime import datetime
# 核心修改：替换为多进程池
from concurrent.futures import ProcessPoolExecutor, as_completed

import logging
# 配置日志：进程安全、实时输出、带时间/进程ID
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(processName)s - %(message)s',  # 改为processName
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.StreamHandler()  # 输出到控制台
    ]
)
logger = logging.getLogger(__name__)


def parse_args():
    """解析命令行参数（仅需输入两个文件夹，无正则）"""
    parser = argparse.ArgumentParser(description='ATAC+甲基化BW全组合合并（笛卡尔积+多进程加速）')
    # 核心输入输出参数（仅需两个文件夹，无需正则）
    parser.add_argument('--atac_dir', '-a', required=True, 
                        help='ATAC BW文件目录（必填，如 ./atac_bw）')
    parser.add_argument('--methyl_dir', '-m', required=True, 
                        help='甲基化 BW文件目录（必填，如 ./methyl_bw）')
    parser.add_argument('--chrom_sizes', '-c', required=True, 
                        help='染色体长度文件路径（必填，如 ./mm10.chrom.sizes）')
    parser.add_argument('--output_dir', '-o', default='./merged_all_pairs', 
                        help='输出目录（默认：./merged_all_pairs）')
    # 计算参数
    parser.add_argument('--bin_size', '-b', type=int, default=10, 
                        help='分辨率（bp，默认10）')
    parser.add_argument('--min_signal', '-t', type=float, default=0.0, 
                        help='最小信号阈值（默认0.0，无过滤）')
    # 并行参数（核心修改：改为进程数）
    parser.add_argument('--n_processes', '-n', type=int, default=4, 
                        help='并行进程数（默认4，建议设为服务器物理核心数）')
    parser.add_argument('--gen_shell', '-g', action='store_true', 
                        help='是否生成批量运行的shell脚本（而非直接运行）')
    parser.add_argument('--shell_out', '-sh', default='run_merge_all_pairs.sh', 
                        help='生成的shell脚本路径（默认：run_merge_all_pairs.sh）')
    
    args = parser.parse_args()
    # 参数合法性检查
    for d in [args.atac_dir, args.methyl_dir]:
        if not os.path.isdir(d):
            logger.info(f"❌ 目录不存在：{d}")
            sys.exit(1)
    if not os.path.isfile(args.chrom_sizes):
        logger.info(f"❌ 染色体文件不存在：{args.chrom_sizes}")
        sys.exit(1)
    return args

def get_all_pairs(atac_dir, methyl_dir):
    """获取ATAC和甲基化文件的所有两两组合（笛卡尔积），可视化输出"""
    # 1. 扫描ATAC文件
    logger.info(f"\n📁 扫描ATAC目录：{atac_dir}")
    atac_files = sorted(glob.glob(os.path.join(atac_dir, "*.bw")))
    atac_names = [os.path.basename(f).replace(".bw", "") for f in atac_files]
    logger.info(f"   找到 {len(atac_files)} 个ATAC BW文件：")
    for i, name in enumerate(atac_names):
        logger.info(f"     [{i+1}] {name}.bw")
    
    # 2. 扫描甲基化文件
    logger.info(f"\n📁 扫描甲基化目录：{methyl_dir}")
    methyl_files = sorted(glob.glob(os.path.join(methyl_dir, "*.bw")))
    methyl_names = [os.path.basename(f).replace(".bw", "") for f in methyl_files]
    logger.info(f"   找到 {len(methyl_files)} 个甲基化 BW文件：")
    for i, name in enumerate(methyl_names):
        logger.info(f"     [{i+1}] {name}.bw")
    
    # 3. 生成所有两两组合（笛卡尔积）
    if len(atac_files) == 0 or len(methyl_files) == 0:
        logger.info("\n❌ ATAC或甲基化目录为空，退出")
        sys.exit(1)
    
    all_pairs = []
    pair_names = []
    for a_idx, (a_file, a_name) in enumerate(zip(atac_files, atac_names)):
        for m_idx, (m_file, m_name) in enumerate(zip(methyl_files, methyl_names)):
            pair_id = f"{a_name}_VS_{m_name}"  # 组合名：ATAC名_VS_甲基化名
            all_pairs.append((pair_id, a_file, m_file))
            pair_names.append(pair_id)
    
    # 4. 输出组合汇总
    total_pairs = len(all_pairs)
    logger.info(f"\n🎯 生成所有两两组合（笛卡尔积）：")
    logger.info(f"   ATAC文件数：{len(atac_files)}")
    logger.info(f"   甲基化文件数：{len(methyl_files)}")
    logger.info(f"   总组合数：{total_pairs} = {len(atac_files)} × {len(methyl_files)}")
    logger.info(f"   前5个组合示例：")
    for i in range(min(5, total_pairs)):
        logger.info(f"     [{i+1}] {pair_names[i]}")
    if total_pairs > 5:
        logger.info(f"     ... 共{total_pairs}个组合")
    
    return all_pairs

def merge_single_pair(params):
    """单个组合合并函数（供多进程调用）"""
    pair_id, atac_bw_path, methyl_bw_path, args = params
    start_time = time.time()
    
    # 创建组合输出目录（避免文件冲突）
    pair_out_dir = os.path.join(args.output_dir, pair_id)
    os.makedirs(pair_out_dir, exist_ok=True)
    
    # 输出文件路径
    bed_path = os.path.join(pair_out_dir, f"{pair_id}_merged.bed")
    log_path = os.path.join(pair_out_dir, f"{pair_id}_merge.log")
    
    # 打开日志文件
    with open(log_path, 'w') as log_f:
        log_f.write(f"开始时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        log_f.write(f"ATAC文件：{os.path.basename(atac_bw_path)}\n")
        log_f.write(f"甲基化文件：{os.path.basename(methyl_bw_path)}\n")
        log_f.write(f"组合名：{pair_id}\n")
        log_f.write(f"分辨率：{args.bin_size}bp\n")
        
        # 打开BW文件
        try:
            atac_bw = pyBigWig.open(atac_bw_path)
            methyl_bw = pyBigWig.open(methyl_bw_path)
            if atac_bw is None or methyl_bw is None:
                err_msg = "无法打开ATAC/甲基化文件"
                log_f.write(f"错误：{err_msg}\n")
                logger.info(f"❌ {pair_id}：{err_msg}")
                return (pair_id, False, err_msg)
        except Exception as e:
            err_msg = f"打开文件失败：{str(e)[:100]}"
            log_f.write(f"错误：{err_msg}\n")
            logger.info(f"❌ {pair_id}：{err_msg}")
            return (pair_id, False, err_msg)
        
        # 读取染色体长度
        chrom_sizes = {}
        with open(args.chrom_sizes, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or '#' in line:
                    continue
                chrom, size = line.split('\t')
                chrom_sizes[chrom] = int(size)
        
        # 写入BED文件（强制保留双信号值，无过滤）
        log_f.write(f"开始写入BED文件：{bed_path}\n")
        processed_bins = 0
        with open(bed_path, 'w') as bed_f:
            for chrom in chrom_sizes.keys():
                # 跳过非标准染色体（减少无效数据）
                if 'random' in chrom or 'chrM' in chrom or 'Un' in chrom:
                    continue
                chrom_size = chrom_sizes[chrom]
                
                # 按bin_size遍历区间
                for start in range(0, chrom_size, args.bin_size):
                    end = min(start + args.bin_size, chrom_size)
                    if start >= end:
                        continue
                    processed_bins += 1
                    
                    # 提取ATAC信号（NaN/Inf替换为0）
                    atac_sig_list = atac_bw.values(chrom, start, end)
                    atac_sig = np.nan if (atac_sig_list is None or len(atac_sig_list) == 0) else np.nanmean(atac_sig_list)
                    atac_sig = 0.0 if np.isnan(atac_sig) or np.isinf(atac_sig) else atac_sig
                    
                    # 提取甲基化信号（NaN/Inf替换为0）
                    methyl_sig_list = methyl_bw.values(chrom, start, end)
                    methyl_sig = np.nan if (methyl_sig_list is None or len(methyl_sig_list) == 0) else np.nanmean(methyl_sig_list)
                    methyl_sig = 0.0 if np.isnan(methyl_sig) or np.isinf(methyl_sig) else methyl_sig
                    
                    # 强制输出5列BED：chr start end ATAC值 甲基化值（无任何过滤）
                    bed_f.write(f"{chrom}\t{start}\t{end}\t{atac_sig:.4f}\t{methyl_sig:.4f}\n")
        
        # 关闭文件
        atac_bw.close()
        methyl_bw.close()
        
        # 计算耗时和文件大小
        cost_time = time.time() - start_time
        bed_size = os.path.getsize(bed_path)/1024/1024 if os.path.exists(bed_path) else 0.0
        
        # 写入日志
        log_f.write(f"完成时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        log_f.write(f"总耗时：{cost_time:.2f}秒\n")
        log_f.write(f"处理区间数：{processed_bins}\n")
        log_f.write(f"输出BED文件大小：{bed_size:.2f}MB\n")
        
        logger.info(f"✅ {pair_id}：完成（耗时{cost_time:.2f}秒 | 处理{processed_bins}个区间）")
        return (pair_id, True, f"完成（耗时{cost_time:.2f}秒 | BED大小{bed_size:.2f}MB）")

def generate_shell_script(all_pairs, args):
    """生成批量运行的shell脚本（笛卡尔积组合）"""
    shell_content = f"""#!/bin/bash
# ATAC+甲基化BW全组合合并脚本（笛卡尔积）
# 生成时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
# 总组合数：{len(all_pairs)}
# 并行建议：根据服务器CPU数调整后台运行数

# 创建输出目录
mkdir -p {args.output_dir}

# 定义单个组合合并函数
merge_pair() {{
    pair_id=$1
    atac_bw=$2
    methyl_bw=$3
    bin_size=$4
    min_signal=$5
    chrom_sizes=$6
    output_dir=$7
    
    # 创建组合输出目录
    mkdir -p $output_dir/$pair_id
    
    # 运行合并（调用当前脚本的单组合模式）
    python {os.path.abspath(__file__)} --single_pair \\
        --pair_id $pair_id \\
        --atac_bw $atac_bw \\
        --methyl_bw $methyl_bw \\
        --chrom_sizes $chrom_sizes \\
        --output_dir $output_dir \\
        --bin_size $bin_size \\
        --min_signal $min_signal > $output_dir/$pair_id/run.log 2>&1
    
    echo "完成：$pair_id"
}}

# 批量运行（每{args.n_processes}个组合后等待，避免过载）
"""
    # 添加所有组合的运行命令
    for i, (pair_id, atac_file, methyl_file) in enumerate(all_pairs):
        # 每N个组合后添加wait（控制并行数）
        if i % args.n_processes == 0 and i > 0:
            shell_content += "wait\n"
        shell_content += f"merge_pair {pair_id} {atac_file} {methyl_file} {args.bin_size} {args.min_signal} {args.chrom_sizes} {args.output_dir} &\n"
    
    # 最后等待所有任务完成
    shell_content += f"""
wait
echo "🎉 所有{len(all_pairs)}个组合合并完成！"
echo "输出目录：{args.output_dir}"
"""
    # 写入shell文件并添加执行权限
    with open(args.shell_out, 'w') as f:
        f.write(shell_content)
    os.chmod(args.shell_out, 0o755)
    
    logger.info(f"\n📜 生成批量shell脚本：{args.shell_out}")
    logger.info(f"   使用方法：")
    logger.info(f"   1. 赋予执行权限：chmod +x {args.shell_out}")
    logger.info(f"   2. 运行脚本：./{args.shell_out}")
    logger.info(f"   3. 并行数控制：脚本默认每{args.n_processes}个组合等待一次，可手动修改")

def main():
    args = parse_args()
    logger.info(f"🚀 开始ATAC+甲基化BW全组合合并（笛卡尔积+多进程）")
    logger.info(f"   运行时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"   ATAC目录：{args.atac_dir}")
    logger.info(f"   甲基化目录：{args.methyl_dir}")
    logger.info(f"   输出目录：{args.output_dir}")
    logger.info(f"   染色体文件：{args.chrom_sizes}")
    logger.info(f"   分辨率：{args.bin_size}bp")
    logger.info(f"   并行进程数：{args.n_processes}")  # 改为进程数
    logger.info(f"   生成shell脚本：{args.gen_shell}")
    
    # 1. 获取所有两两组合（笛卡尔积）
    all_pairs = get_all_pairs(args.atac_dir, args.methyl_dir)
    total_pairs = len(all_pairs)
    if total_pairs == 0:
        logger.info("\n❌ 无可用组合，退出")
        sys.exit(1)
    
    # 2. 生成shell脚本（或直接多进程运行）
    if args.gen_shell:
        generate_shell_script(all_pairs, args)
        return
    
    # 3. 多进程并行运行所有组合（核心修改）
    logger.info(f"\n⚡ 开始多进程合并（进程数：{args.n_processes} | 总组合数：{total_pairs}）")
    start_total = time.time()
    
    # 准备多进程参数
    params_list = [(pair_id, a_file, m_file, args) for pair_id, a_file, m_file in all_pairs]
    
    # 启动进程池（核心修改：替换为ProcessPoolExecutor）
    with ProcessPoolExecutor(max_workers=args.n_processes) as executor:
        future_to_pair = {executor.submit(merge_single_pair, p): p[0] for p in params_list}
        results = []
        # 实时收集结果
        for future in as_completed(future_to_pair):
            pair_id = future_to_pair[future]
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                err_msg = f"运行出错：{str(e)[:100]}"
                results.append((pair_id, False, err_msg))
                logger.info(f"❌ {pair_id}：{err_msg}")
    
    # 输出最终汇总
    total_time = time.time() - start_total
    success = [r for r in results if r[1]]
    fail = [r for r in results if not r[1]]
    
    logger.info(f"\n📊 最终汇总（总耗时：{total_time:.2f}秒）：")
    logger.info(f"   总组合数：{total_pairs}")
    logger.info(f"   成功数：{len(success)}")
    logger.info(f"   失败数：{len(fail)}")
    logger.info(f"   平均耗时：{total_time/total_pairs:.2f}秒/组合")
    if fail:
        logger.info(f"   失败组合：{[r[0] for r in fail]}")
    logger.info(f"\n🎉 所有任务完成！输出目录：{args.output_dir}")

if __name__ == "__main__":
    # 支持单组合运行模式（供shell脚本调用）
    if '--single_pair' in sys.argv:
        single_parser = argparse.ArgumentParser()
        single_parser.add_argument('--single_pair', action='store_true')
        single_parser.add_argument('--pair_id', required=True)
        single_parser.add_argument('--atac_bw', required=True)
        single_parser.add_argument('--methyl_bw', required=True)
        single_parser.add_argument('--chrom_sizes', required=True)
        single_parser.add_argument('--output_dir', required=True)
        single_parser.add_argument('--bin_size', type=int, default=200)
        single_parser.add_argument('--min_signal', type=float, default=0.0)
        single_args = single_parser.parse_args()
        merge_single_pair((single_args.pair_id, single_args.atac_bw, single_args.methyl_bw, single_args))
    else:
        main()