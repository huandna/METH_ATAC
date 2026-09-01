import os
import re
import glob
import numpy as np
import subprocess
import pyBigWig  # 核心依赖：专门处理BigWig文件
from collections import defaultdict

def get_chrom_sizes(chrom_sizes_path):
    """读取染色体长度文件，返回字典{染色体: 长度}"""
    chrom_sizes = {}
    with open(chrom_sizes_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or '#' in line:
                continue
            chrom, size = line.split('\t')
            chrom_sizes[chrom] = int(size)
    return chrom_sizes

def get_celltype_groups(bw_dir):
    """按文件名分组：适配 celltype_数字.bw 格式（如 Neural_12.bw）"""
    bw_files = glob.glob(os.path.join(bw_dir, "*.bw"))
    celltype_groups = defaultdict(list)
    
    # 修正后的正则：匹配 celltype_数字.bw
    pattern = re.compile(r'^(.*?)_\d+\.bw$')
    for bw_file in bw_files:
        fname = os.path.basename(bw_file)
        match = pattern.match(fname)
        if match:
            celltype = match.group(1)
            celltype_groups[celltype].append(bw_file)
    
    # 检查分组结果
    print("=== 细胞类型分组结果 ===")
    if not celltype_groups:
        print("⚠️ 未找到符合格式的BW文件（需为 celltype_数字.bw）")
    else:
        for celltype, files in celltype_groups.items():
            print(f"{celltype}: {len(files)}个重复文件")
            for f in files:
                print(f"  - {os.path.basename(f)}")
    print("========================")
    return celltype_groups

def merge_bw_files(celltype, bw_files, chrom_sizes, output_dir, bin_size=10):
    """整合单个细胞类型的所有BW文件，输出bedGraph和BW"""
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 用pyBigWig打开所有BW文件
    bw_handles = []
    for f in bw_files:
        try:
            bw = pyBigWig.open(f)
            bw_handles.append(bw)
        except Exception as e:
            print(f"⚠️ 打开文件 {f} 失败：{e}")
            continue
    
    if len(bw_handles) == 0:
        print(f"❌ {celltype} 无可用的BW文件，跳过")
        return
    
    # 输出bedGraph文件路径
    bedgraph_path = os.path.join(output_dir, f"{celltype}.bedGraph")
    
    # 逐染色体处理
    with open(bedgraph_path, 'w') as bg_f:
        for chrom in chrom_sizes.keys():
            # 跳过非标准染色体（减少计算量）
            if 'random' in chrom or 'chrM' in chrom or 'Un' in chrom:
                continue
            chrom_size = chrom_sizes[chrom]
            print(f"处理 {celltype} - {chrom} (长度: {chrom_size}bp)")
            
            # 按bin_size分块遍历（避免内存溢出）
            for start in range(0, chrom_size, bin_size):
                end = min(start + bin_size, chrom_size)
                if start >= end:
                    continue
                
                # 提取所有重复文件在该区间的信号
                signals = []
                for bw in bw_handles:
                    try:
                        # pyBigWig提取信号（返回list，需转numpy）
                        sig_list = bw.values(chrom, start, end)
                        if sig_list is None:
                            continue
                        # 转为numpy数组并过滤NaN/Inf
                        sig = np.array(sig_list, dtype=np.float64)
                        sig = sig[~np.isnan(sig)]
                        sig = sig[~np.isinf(sig)]
                        if len(sig) == 0:
                            continue
                        # 计算该bin的均值
                        bin_mean = np.mean(sig)
                        if bin_mean > 0:  # 过滤无信号位点
                            signals.append(bin_mean)
                    except (KeyError, ValueError, TypeError):
                        # 跳过无信号/不存在的区间
                        continue
                
                # 计算所有重复的均值并写入bedGraph
                if len(signals) > 0:
                    final_mean = np.mean(signals)
                    # 保留4位小数，符合bedGraph格式规范
                    bg_f.write(f"{chrom}\t{start}\t{end}\t{final_mean:.4f}\n")
    
    # 关闭所有BW句柄
    for bw in bw_handles:
        bw.close()
    
    # 将bedGraph转为BW文件（调用UCSC工具）
    bw_path = os.path.join(output_dir, f"{celltype}.bw")
    chrom_sizes_path = os.path.join(output_dir, f"{celltype}_chrom.sizes")
    
    # 保存过滤后的染色体长度文件（仅保留有效染色体）
    with open(chrom_sizes_path, 'w') as f:
        for chrom, size in chrom_sizes.items():
            if 'random' not in chrom and 'chrM' not in chrom and 'Un' not in chrom:
                f.write(f"{chrom}\t{size}\n")
    
    # 执行bedGraphToBigWig命令
    cmd = [
        "bedGraphToBigWig",
        bedgraph_path,
        chrom_sizes_path,
        bw_path
    ]
    try:
        # 执行命令并捕获输出
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        print(f"✅ {celltype} 整合完成：")
        print(f"  - bedGraph: {bedgraph_path}")
        print(f"  - BW: {bw_path}")
    except subprocess.CalledProcessError as e:
        print(f"❌ {celltype} BW转换失败：")
        print(f"  错误信息：{e.stderr}")
    except FileNotFoundError:
        print(f"❌ 未找到bedGraphToBigWig工具，请安装：")
        print(f"     conda install -c bioconda ucsc-bedgraphtobigwig -y")
    finally:
        # 清理临时的染色体长度文件
        if os.path.exists(chrom_sizes_path):
            os.remove(chrom_sizes_path)

def main():
    """主函数：配置参数并执行整合"""
    # ====================== 请手动配置以下参数 ======================
    BW_DIR = "./bw"          # BW文件所在目录（修改为你的实际路径）
    OUTPUT_DIR = "./"# 输出目录（自动创建）
    CHROM_SIZES_PATH = "./species_chrom_size.xls"  # 染色体长度文件（修改为你的版本）
    BIN_SIZE = 100                  # 分辨率（建议10/50，减少文件体积）
    # ===============================================================
    
    # 1. 检查pyBigWig是否安装
    try:
        import pyBigWig
    except ImportError:
        print("❌ 未安装pyBigWig，正在自动安装...")
        subprocess.run(["pip", "install", "pyBigWig", "-i", "https://pypi.tuna.tsinghua.edu.cn/simple"], check=True)
        import pyBigWig
    
    # 2. 读取染色体长度文件
    print("读取染色体长度文件...")
    try:
        chrom_sizes = get_chrom_sizes(CHROM_SIZES_PATH)
        print(f"✅ 成功读取 {len(chrom_sizes)} 条染色体信息")
    except FileNotFoundError:
        print(f"❌ 染色体长度文件不存在：{CHROM_SIZES_PATH}")
        print("请下载对应基因组版本的chrom.sizes文件：")
        print("  mm10: wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes")
        print("  hg38: wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes")
        return
    
    # 3. 按细胞类型分组BW文件
    print("\n分组BW文件...")
    celltype_groups = get_celltype_groups(BW_DIR)
    if not celltype_groups:
        return
    
    # 4. 逐个整合细胞类型
    print("\n开始整合...")
    for celltype, bw_files in celltype_groups.items():
        merge_bw_files(
            celltype=celltype,
            bw_files=bw_files,
            chrom_sizes=chrom_sizes,
            output_dir=OUTPUT_DIR,
            bin_size=BIN_SIZE
        )
    
    print("\n🎉 所有细胞类型整合完成！输出目录：", OUTPUT_DIR)

if __name__ == "__main__":
    main()