import os
import gzip
import pandas as pd
import re
import glob
import matplotlib.pyplot as plt
import numpy as np
def process_methylation_files(file_pattern):
    """
    处理甲基化数据文件，统计每个样本在不同上下文类型下的甲基化程度分布
    
    参数:
        file_pattern: 文件匹配模式
    
    返回:
        包含甲基化程度分布统计的DataFrame
    """
    files = glob.glob(file_pattern)
    print(f"找到 {len(files)} 个匹配的文件")
    
    if not files:
        print("警告：没有找到匹配的文件")
        return None
    
    # 定义甲基化值区间 (0-1, 每0.1一个区间)
    bins = np.arange(0, 1.1, 0.1)
    bin_labels = [f"{i*10}" for i in range(len(bins)-1)]
    
    all_data = []
    
    for i, file in enumerate(files):
        print(f"\n处理文件 {i+1}/{len(files)}: {file}")
        
        filename = os.path.basename(file)
        sample_name = filename.replace('.mC_level_Identification_stat.txt.gz', '')
        print(f"样本名: {sample_name}")
        
        try:
            with gzip.open(file, 'rt') as f:
                df = pd.read_csv(f, sep='\t', header=None, 
                               names=['Chromosome', 'Position', 'Strand', 'mC_counts', 'umC_counts', 
                                      'methylation_level', 'Context', 'Pvalue', 'Corrected_pvalue'])
                
                print(f"读取到 {len(df)} 行数据")
            
            # 按上下文类型分组统计
            for context_type in ['CG', 'CHG', 'CHH']:
                # 筛选当前上下文类型的数据
                context_df = df[df['Context'] == context_type]
                
                if len(context_df) == 0:
                    print(f"样本 {sample_name} 中没有 {context_type} 类型的数据")
                    continue
                
                # 统计甲基化值分布
                methylation_levels = context_df['methylation_level']
                counts, _ = np.histogram(methylation_levels, bins=bins)
                
                # 计算频率
                frequencies = counts / counts.sum() * 100
                
                # 创建结果DataFrame
                sample_data = pd.DataFrame({
                    'Sample': sample_name,
                    'Context_Type': context_type,
                    'Methylation_Level': bin_labels,
                    'Frequency': frequencies
                })
                
                all_data.append(sample_data)
            
            print(f"样本 {sample_name} 处理完成")
            
        except Exception as e:
            print(f"处理文件 {file} 时出错: {str(e)}")
            continue
    
    if not all_data:
        print("警告：没有成功处理任何文件")
        return None
    
    # 合并所有样本数据
    result_df = pd.concat(all_data, ignore_index=True)
    print(f"\n合并后的数据形状: {result_df.shape}")
    
    return result_df

def process_replicates(df):
    """
    处理重复样本，取平均值
    
    参数:
        df: 包含所有样本数据的DataFrame
    
    返回:
        处理重复样本后的DataFrame
    """
    # 提取样本基础名称（去除末尾的数字）
    df['Base_Sample'] = df.index.str.replace(r'\d+$', '', regex=True)
    
    # 按基础样本名分组，计算平均值
    averaged_df = df.groupby('Base_Sample').mean()
    
    return averaged_df
if __name__ == "__main__":
    # 替换为你的文件路径模式
    file_pattern = "/home/yaoxinw/workdir/data/04.chidian_5mc/Result_X101SC25117400-Z01-J001_Epinephelus_akaara/Result_X101SC25117400-Z01-J001_Epinephelus_akaara/3.MethylationStat/3.1.MLevel_stat/*.mC_level_Identification_stat.txt.gz"
    
    # 处理数据
    methylation_distribution = process_methylation_files(file_pattern)
    
    if methylation_distribution is not None:
        # 保存原始数据
        methylation_distribution.to_csv("methylation_distribution_raw.csv", index=False)
        print("\n数据处理完成，原始数据已保存到 methylation_distribution_raw.csv")
        
        # 保存合并数据
        pivot_df = methylation_distribution.pivot_table(
            index=['Methylation_Level', 'Context_Type'], 
            columns='Sample', 
            values='Frequency'
        ).reset_index()
        
        pivot_df.to_csv('methylation_distribution_all.csv', index=False)
        print("合并的甲基化分布数据已保存到 methylation_distribution_all.csv")
    else:
        print("\n数据处理失败，请检查错误信息")

