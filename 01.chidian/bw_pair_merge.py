import os
import glob
import pyBigWig
import numpy as np

# ====================== 你只需要改这里 ======================
atac_dir = "./atac"        # ATAC的bw文件夹
meth_dir = "./methyl"      # 甲基化bw文件夹
out_dir  = "./merged_bed"  # 输出目录
bin_size = 10              # 窗口大小
chrom_file = "./mm10.chrom.sizes"  # 染色体长度文件
# ==========================================================

os.makedirs(out_dir, exist_ok=True)

# 读取染色体
chroms = {}
with open(chrom_file) as f:
    for line in f:
        c, s = line.strip().split()
        chroms[c] = int(s)

# 匹配样本（按文件名前缀配对）
atac_dict = {os.path.basename(f).rsplit('_',1)[0]: f for f in glob.glob(f"{atac_dir}/*.bw")}
meth_dict = {os.path.basename(f).rsplit('_',1)[0]: f for f in glob.glob(f"{meth_dir}/*.bw")}

samples = set(atac_dict.keys()) & set(meth_dict.keys())

if not samples:
    print("❌ 没有匹配到成对样本")
    exit()

print(f"✅ 匹配到 {len(samples)} 个样本")

# 逐个合并
for sample in sorted(samples):
    atac_bw = pyBigWig.open(atac_dict[sample])
    meth_bw = pyBigWig.open(meth_dict[sample])

    out_bed = os.path.join(out_dir, f"{sample}.bed")
    print(f"正在处理：{sample}")

    with open(out_bed, 'w') as bed:
        for chrom, size in chroms.items():
            if '_' in chrom or 'chrM' in chrom:
                continue

            for start in range(0, size, bin_size):
                end = start + bin_size
                if end > size:
                    end = size

                # 取信号
                a = atac_bw.values(chrom, start, end)
                m = meth_bw.values(chrom, start, end)

                a_val = float(np.nanmean(a)) if a else 0.0
                m_val = float(np.nanmean(m)) if m else 0.0

                # 直接输出四列：chr start end 信号
                # bed.write(f"{chrom}\t{start}\t{end}\t{a_val:.4f}\n")
                # 五列：chr start end ATAC 甲基化
                bed.write(f"{chrom}\t{start}\t{end}\t{a_val:.4f}\t{m_val:.4f}\n")

    atac_bw.close()
    meth_bw.close()

print("\n🎉 全部完成！信号列一定存在")