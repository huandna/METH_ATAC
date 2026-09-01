import os
import re
import glob
import logging
import numpy as np
import pyBigWig
import subprocess
import argparse
from collections import defaultdict
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
# 替换原有的 ThreadPoolExecutor 导入，新增 ProcessPoolExecutor
from concurrent.futures import ProcessPoolExecutor, as_completed

# ====================== 日志配置（线程安全）======================
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(threadName)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.StreamHandler(),  # 输出到控制台
        logging.FileHandler('merge_methyl_replicates.log', encoding='utf-8')  # 输出到日志文件
    ]
)
logger = logging.getLogger(__name__)

def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='甲基化BW重复样本合并（多线程+日志版）')
    parser.add_argument('--input_dir', '-i', required=True,
                        help='输入BW文件目录（包含所有重复样本，如 ./methyl_bw）')
    parser.add_argument('--output_dir', '-o', default='./merged_methyl_replicates',
                        help='输出目录（默认：./merged_methyl_replicates）')
    parser.add_argument('--chrom_sizes', '-c', required=True,
                        help='染色体长度文件（如 ./mm10.chrom.sizes）')
    parser.add_argument('--bin_size', '-b', type=int, default=1,
                        help='分辨率（bp，甲基化建议1bp）')
    parser.add_argument('--min_signal', '-t', type=float, default=0.0,
                        help='最小信号阈值（默认0，无过滤）')
    parser.add_argument('--n_threads', '-n', type=int, default=4,
                        help='并行线程数（默认4，建议设为CPU核心数的50%-80%）')
    return parser.parse_args()

def group_replicates(input_dir):
    """按「组织名+修饰类型」分组重复样本（返回有效分组）"""
    #pattern = re.compile(r'^([a-zA-Z]+)\d+\.([CAG]+)\.bw$')
    pattern = re.compile(r'^([a-zA-Z]+)\d+\.(CG|CHG|CHH)\.bw$')
    file_groups = defaultdict(list)
    
    # 扫描所有BW文件
    bw_files = glob.glob(os.path.join(input_dir, "*.bw"))
    logger.info(f"扫描到 {len(bw_files)} 个BW文件")
    
    for f in bw_files:
        fname = os.path.basename(f)
        match = pattern.match(fname)
        if match:
            tissue = match.group(1)
            mod_type = match.group(2)
            group_name = f"{tissue}.{mod_type}"
            file_groups[group_name].append(f)
            logger.debug(f"匹配文件：{fname} → 分组 {group_name}")
        else:
            logger.warning(f"文件不匹配规则，跳过：{fname}")
    
    # 验证分组（每组需3个重复）
    valid_groups = {}
    for group_name, files in file_groups.items():
        if len(files) == 3:
            valid_groups[group_name] = files
            logger.info(f"有效分组：{group_name} → {len(files)}个重复")
            for f in files:
                logger.debug(f"  - {os.path.basename(f)}")
        else:
            logger.error(f"无效分组：{group_name} → 仅{len(files)}个重复（需3个）")
    
    # 分组汇总
    logger.info(f"分组汇总：总分组数{len(file_groups)} | 有效分组{len(valid_groups)} | 无效分组{len(file_groups)-len(valid_groups)}")
    
    if not valid_groups:
        logger.critical("无有效重复分组，退出脚本")
        exit(1)
    return valid_groups

def merge_group(params):
    """
    单个分组合并函数（供多线程调用）
    :param params: (group_name, bw_files, chrom_sizes_path, output_dir, bin_size, min_signal)
    """
    group_name, bw_files, chrom_sizes_path, output_dir, bin_size, min_signal = params
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 输出文件路径
    merged_bed_path = os.path.join(output_dir, f"{group_name}.bed")
    merged_bg_path = os.path.join(output_dir, f"{group_name}.bedGraph")
    merged_bw_path = os.path.join(output_dir, f"{group_name}.bw")
    
    logger.info(f"开始合并分组：{group_name}")
    
    try:
        # 打开所有重复BW文件
        bw_handles = []
        for f in bw_files:
            bw = pyBigWig.open(f)
            if bw is None:
                logger.error(f"无法打开文件：{os.path.basename(f)}")
                raise Exception(f"文件打开失败：{os.path.basename(f)}")
            bw_handles.append(bw)
        
        # 读取染色体长度
        chrom_sizes = {}
        with open(chrom_sizes_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or '#' in line:
                    continue
                chrom, size = line.split('\t')
                chrom_sizes[chrom] = int(size)
        
        # 合并信号（取3个重复的均值）
        processed_chroms = 0
        processed_bins = 0
        with open(merged_bed_path, 'w') as bed_f, open(merged_bg_path, 'w') as bg_f:
            for chrom in chrom_sizes.keys():
                if 'random' in chrom or 'chrM' in chrom or 'Un' in chrom:
                    continue
                chrom_size = chrom_sizes[chrom]
                processed_chroms += 1
                
                for start in range(0, chrom_size, bin_size):
                    end = min(start + bin_size, chrom_size)
                    if start >= end:
                        continue
                    processed_bins += 1
                    
                    # 提取3个重复的信号值
                    replicate_signals = []
                    for bw in bw_handles:
                        sig_list = bw.values(chrom, start, end)
                        if sig_list is None or len(sig_list) == 0:
                            continue
                        sig_array = np.array(sig_list, dtype=np.float64)
                        sig_array = np.nan_to_num(sig_array, nan=0.0, posinf=0.0, neginf=0.0)
                        rep_mean = np.mean(sig_array)
                        replicate_signals.append(rep_mean)
                    
                    # 计算均值并过滤
                    if len(replicate_signals) == 3:
                        merged_mean = np.mean(replicate_signals)
                        if merged_mean >= min_signal:
                            bed_f.write(f"{chrom}\t{start}\t{end}\t{merged_mean:.4f}\n")
                            bg_f.write(f"{chrom}\t{start}\t{end}\t{merged_mean:.4f}\n")
        
        # 关闭BW句柄
        for bw in bw_handles:
            bw.close()
        
        # bedGraph转BW
        cmd = ["bedGraphToBigWig", merged_bg_path, chrom_sizes_path, merged_bw_path]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # 清理临时文件
        if os.path.exists(merged_bg_path):
            os.remove(merged_bg_path)
        
        logger.info(f"分组 {group_name} 合并成功：")
        logger.info(f"  - BED文件：{os.path.basename(merged_bed_path)}")
        logger.info(f"  - BW文件：{os.path.basename(merged_bw_path)}")
        logger.info(f"  - 处理染色体数：{processed_chroms} | 处理区间数：{processed_bins}")
        
        return (group_name, True)
    
    except subprocess.CalledProcessError as e:
        logger.error(f"分组 {group_name} BW转换失败：{e.stderr[:100]}")
        return (group_name, False)
    except FileNotFoundError:
        logger.error(f"分组 {group_name} 缺少bedGraphToBigWig工具，请安装：conda install -c bioconda ucsc-bedgraphtobigwig")
        return (group_name, False)
    except Exception as e:
        logger.error(f"分组 {group_name} 合并失败：{str(e)[:100]}", exc_info=True)
        return (group_name, False)

def main():
    args = parse_args()
    logger.info("="*50)
    logger.info("甲基化重复样本合并脚本（多线程+日志版）启动")
    logger.info(f"启动时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"输入目录：{args.input_dir}")
    logger.info(f"输出目录：{args.output_dir}")
    logger.info(f"染色体文件：{args.chrom_sizes}")
    logger.info(f"分辨率：{args.bin_size}bp | 最小信号阈值：{args.min_signal}")
    logger.info(f"并行线程数：{args.n_threads}")
    logger.info("="*50)
    
    # 1. 分组重复样本
    valid_groups = group_replicates(args.input_dir)
    
    # 2. 准备多线程参数
    params_list = []
    for group_name, files in valid_groups.items():
        params = (
            group_name,
            files,
            args.chrom_sizes,
            args.output_dir,
            args.bin_size,
            args.min_signal
        )
        params_list.append(params)
    
    # 3. 多线程并行合并
    logger.info(f"启动多线程合并，共 {len(params_list)} 个分组，线程数 {args.n_threads}")
    success_count = 0
    fail_count = 0
    failed_groups = []
    # 改后多进程代码（绕过GIL，利用多核）
    with ProcessPoolExecutor(max_workers=args.n_threads) as executor:
        future_to_group = {executor.submit(merge_group, p): p[0] for p in params_list}
        # 提交任务
        
        # 收集结果
        for future in as_completed(future_to_group):
            group_name = future_to_group[future]
            try:
                group_name, result = future.result()
                if result:
                    success_count += 1
                else:
                    fail_count += 1
                    failed_groups.append(group_name)
            except Exception as e:
                logger.error(f"获取分组 {group_name} 结果失败：{str(e)[:50]}")
                fail_count += 1
                failed_groups.append(group_name)
    
    # 4. 最终汇总
    logger.info("="*50)
    logger.info("合并任务完成汇总：")
    logger.info(f"总分组数：{len(valid_groups)}")
    logger.info(f"成功合并：{success_count}")
    logger.info(f"合并失败：{fail_count}")
    if failed_groups:
        logger.error(f"失败分组列表：{failed_groups}")
    logger.info(f"输出目录：{args.output_dir}")
    logger.info(f"日志文件：merge_methyl_replicates.log")
    logger.info("="*50)

if __name__ == "__main__":
    main()