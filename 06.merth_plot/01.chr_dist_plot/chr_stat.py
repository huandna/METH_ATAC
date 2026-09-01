#!/usr/bin/env python3
import os
import gzip
import pandas as pd
import re
import glob

def process_methylation_files(file_pattern):
    """
    处理甲基化数据文件，统计每个样本在每条染色体上的甲基化水平
    
    参数:
        file_pattern: 文件匹配模式，如 "./3.1.MLevel_stat/*_level_Identification_stat.txt.gz"
    
    返回:
        包含统计结果的DataFrame
    """
    # 获取所有匹配的文件
    files = glob.glob(file_pattern)
    print(f"找到 {len(files)} 个匹配的文件")
    
    if not files:
        print("警告：没有找到匹配的文件，请检查文件路径和模式是否正确")
        return None
    
    # 存储所有样本的数据
    all_data = []
    
    for i, file in enumerate(files):
        print(f"\n处理文件 {i+1}/{len(files)}: {file}")
        
        # 从完整路径中提取文件名
        filename = os.path.basename(file)
        print(f"文件名: {filename}")
        
        # 提取样本名：去除.mC_level_Identification_stat.txt.gz后缀
        sample_name = filename.replace('.mC_level_Identification_stat.txt.gz', '')
        print(f"样本名: {sample_name}")
        
        try:
            # 读取gzip压缩文件
            with gzip.open(file, 'rt') as f:
                # 读取数据
                df = pd.read_csv(f, sep='\t', header=None, 
                                 names=['Chromosome', 'Position', 'Strand', 'mC_counts', 'umC_counts', 
                                        'methylation_level', 'Context', 'Pvalue', 'Corrected_pvalue'])
                
                print(f"读取到 {len(df)} 行数据")
            
            # 按染色体分组，计算甲基化水平平均值（不筛选染色体）
            chr_stats = df.groupby('Chromosome')['methylation_level'].mean().reset_index()
            
            # 添加样本名
            chr_stats['Sample'] = sample_name
            
            # 添加到结果列表
            all_data.append(chr_stats)
            print(f"样本 {sample_name} 处理完成，统计了 {len(chr_stats)} 条染色体")
            
        except Exception as e:
            print(f"处理文件 {file} 时出错: {str(e)}")
            continue
    
    if not all_data:
        print("警告：没有成功处理任何文件")
        return None
    
    # 合并所有样本数据
    result_df = pd.concat(all_data)
    print(f"\n合并后的数据形状: {result_df.shape}")
    
    # 透视表：行=样本，列=染色体，值=甲基化水平
    pivot_df = result_df.pivot(index='Sample', columns='Chromosome', values='methylation_level')
    
    # 筛选常染色体 (Chr加数字)
    # 使用更灵活的正则表达式，匹配Chr后跟一个或多个数字的格式
    chr_cols = [col for col in pivot_df.columns if re.match(r'^Chr\d+$', col)]
    pivot_df = pivot_df[chr_cols]
    
    print(f"筛选后透视表形状: {pivot_df.shape}")
    print(f"透视表行名: {pivot_df.index.tolist()}")
    print(f"透视表列名: {pivot_df.columns.tolist()}")
    
    return pivot_df

# 使用示例
if __name__ == "__main__":
    # 替换为你的文件路径模式
    file_pattern = "/home/yaoxinw/workdir/data/04.chidian_5mc/Result_X101SC25117400-Z01-J001_Epinephelus_akaara/Result_X101SC25117400-Z01-J001_Epinephelus_akaara/3.MethylationStat/3.1.MLevel_stat/*.mC_level_Identification_stat.txt.gz"
    
    # 处理数据
    methylation_stats = process_methylation_files(file_pattern)
    
    if methylation_stats is not None:
        # 保存结果
        methylation_stats.to_csv("methylation_level_by_chromosome.csv")
        print("\n数据处理完成，结果已保存到 methylation_level_by_chromosome.csv")
    else:
        print("\n数据处理失败，请检查错误信息")
