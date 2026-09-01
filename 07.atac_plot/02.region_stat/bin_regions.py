import pandas as pd
import argparse
import numpy as np

def bin_regions(input_file, chrom_size_file, bin_size=500):
    """将染色体区域按固定bin大小划分并计算信号均值"""
    df = pd.read_csv(input_file)
    
    # 读取染色体长度文件
    chrom_sizes = pd.read_csv(chrom_size_file, sep='\t', header=None, names=['chrom', 'size'])
    chrom_size_dict = dict(zip(chrom_sizes['chrom'], chrom_sizes['size']))
    
    results = []
    
    # 按染色体分组
    for chrom, chrom_group in df.groupby('chrom'):
        # 获取该染色体的实际长度
        chrom_length = chrom_size_dict.get(chrom, chrom_group['end'].max())
        
        # 计算需要的bin数量
        num_bins = max(1, int(np.ceil(chrom_length / bin_size)))
        
        # 为每个bin计算均值
        for bin_idx in range(num_bins):
            bin_start = bin_idx * bin_size + 1  # 从1号碱基开始
            bin_end = (bin_idx + 1) * bin_size
            
            # 创建bin编号，格式为Chr01_1
            bin_id = f"{chrom}_{bin_idx + 1}"
            
            # 找到与当前bin重叠的区域
            mask = ((chrom_group['start'] < bin_end) & (chrom_group['end'] > bin_start))
            overlapping = chrom_group[mask]
            
            if len(overlapping) > 0:
                # 计算加权平均
                weights = []
                values = []
                gene_names = []
                features = []
                
                for _, row in overlapping.iterrows():
                    overlap_start = max(row['start'], bin_start)
                    overlap_end = min(row['end'], bin_end)
                    overlap_length = overlap_end - overlap_start
                    
                    weights.append(overlap_length)
                    values.append(row['signal'])
                    gene_names.append(row['gene_name'])
                    features.append(row['feature'])
                
                avg_signal = np.average(values, weights=weights)
                
                # 记录该bin中最常见的gene_name和feature
                most_common_gene = max(set(gene_names), key=gene_names.count) if gene_names else 'NA'
                most_common_feature = max(set(features), key=features.count) if features else 'NA'
                
                results.append({
                    'bin_id': bin_id,
                    'signal': avg_signal,
                    'gene_name': most_common_gene,
                    'feature': most_common_feature
                })
            else:
                # 如果该bin没有重叠区域，仍然记录该bin，信号值为0
                results.append({
                    'bin_id': bin_id,
                    'signal': 0,
                    'gene_name': 'NA',
                    'feature': 'NA'
                })
    
    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser(description='Bin chromosomal regions and calculate signal averages')
    parser.add_argument('--input', required=True, help='Input CSV file')
    parser.add_argument('--output', required=True, help='Output CSV file')
    parser.add_argument('--bin', type=int, default=500, help='Bin size in bp (default: 500)')
    parser.add_argument('--chrom_size', required=True, help='Chromosome size file (tab-separated)')
    
    args = parser.parse_args()
    
    # 处理binning
    binned_df = bin_regions(args.input, args.chrom_size, args.bin)
    
    # 保存结果
    binned_df.to_csv(args.output, index=False)
    print(f"Binned data saved to {args.output}")

if __name__ == '__main__':
    main()
