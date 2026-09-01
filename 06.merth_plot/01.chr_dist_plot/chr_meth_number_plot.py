#!/usr/bin/env python3
import os
import gzip
import pandas as pd
import re
import glob
import matplotlib.pyplot as plt
import numpy as np

def process_methylation_files(file_pattern):
    """
    处理甲基化数据文件，统计每个样本在每条染色体上甲基化值超过0.5的位点个数
    
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
            
            # 筛选甲基化值超过0.5的位点
            df_filtered = df[df['methylation_level'] > 0.5]
            print(f"甲基化值超过0.5的位点数: {len(df_filtered)}")
            
            # 按染色体分组，计算甲基化值超过0.5的位点个数
            chr_stats = df_filtered.groupby('Chromosome').size().reset_index(name='count')
            
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
    
    # 透视表：行=样本，列=染色体，值=甲基化值超过0.5的位点个数
    pivot_df = result_df.pivot(index='Sample', columns='Chromosome', values='count')
    
    # 筛选常染色体 (Chr加数字)
    chr_cols = [col for col in pivot_df.columns if re.match(r'^Chr\d+$', col)]
    pivot_df = pivot_df[chr_cols]
    
    print(f"筛选后透视表形状: {pivot_df.shape}")
    
    return pivot_df

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

def plot_methylation_stats(df, output_file='methylation_stats.png'):
    """
    绘制甲基化统计折线图
    
    参数:
        df: 包含统计数据的DataFrame
        output_file: 输出图片文件名
    """
    # 转置数据，使染色体为行，样本为列
    df_t = df.T
    
    # 按染色体编号排序
    df_t = df_t.reindex(sorted(df_t.index, key=lambda x: int(x.replace('Chr', ''))))
    
    # 创建折线图
    plt.figure(figsize=(12, 6))
    
    # 为每个样本绘制一条折线
    for sample in df_t.columns:
        plt.plot(df_t.index, df_t[sample], marker='o', label=sample)
    
    # 设置图表标题和标签
    plt.title('Number of Methylation Sites (>0.5) by Chromosome')
    plt.xlabel('Chromosome')
    plt.ylabel('Count of Methylation Sites')
    
    # 添加图例
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"折线图已保存到 {output_file}")
    
    # 显示图表
    plt.show()

# 使用示例
if __name__ == "__main__":
    # 替换为你的文件路径模式
    file_pattern = "/home/yaoxinw/workdir/data/04.chidian_5mc/Result_X101SC25117400-Z01-J001_Epinephelus_akaara/Result_X101SC25117400-Z01-J001_Epinephelus_akaara/3.MethylationStat/3.1.MLevel_stat/*.mC_level_Identification_stat.txt.gz"
    
    # 处理数据
    methylation_stats = process_methylation_files(file_pattern)
    
    if methylation_stats is not None:
        # 保存结果
        methylation_stats.to_csv("methylation_count_by_chromosome.csv")
        print("\n数据处理完成，结果已保存到 methylation_count_by_chromosome.csv")
        
        # 处理重复样本
        averaged_stats = process_replicates(methylation_stats)
        averaged_stats.to_csv("methylation_count_by_chromosome_averaged.csv")
        print("重复样本处理完成，结果已保存到 methylation_count_by_chromosome_averaged.csv")
        
        # 绘制折线图
        plot_methylation_stats(averaged_stats)
    else:
        print("\n数据处理失败，请检查错误信息")