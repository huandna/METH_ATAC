import pandas as pd
import argparse
import os
from collections import defaultdict

def parse_gff(gff_file):
    """解析GFF文件，提取基因信息并建立层级关系"""
    genes = {}
    mrnas = defaultdict(list)
    exons = defaultdict(list)
    cdss = defaultdict(list)
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9:
                chrom, feature, start, end = parts[0], parts[2], int(parts[3]), int(parts[4])
                attributes = parts[8]
                
                # 提取ID和Parent信息
                feature_id = None
                parent_id = None
                for attr in attributes.split(';'):
                    if attr.startswith('ID='):
                        feature_id = attr.split('=')[1]
                    elif attr.startswith('Parent='):
                        parent_id = attr.split('=')[1]
                
                if feature == 'gene' and feature_id:
                    genes[feature_id] = {'chrom': chrom, 'start': start, 'end': end}
                elif feature == 'mRNA' and parent_id:
                    mrnas[parent_id].append({'chrom': chrom, 'start': start, 'end': end, 'id': feature_id})
                elif feature == 'exon' and parent_id:
                    exons[parent_id].append({'chrom': chrom, 'start': start, 'end': end})
                elif feature == 'CDS' and parent_id:
                    cdss[parent_id].append({'chrom': chrom, 'start': start, 'end': end})
    
    # 构建完整的基因结构
    gene_dict = defaultdict(list)
    
    for gene_id, gene_info in genes.items():
        chrom = gene_info['chrom']
        gene_start = gene_info['start']
        gene_end = gene_info['end']
        
        # 添加gene特征
        gene_dict[chrom].append({
            'start': gene_start,
            'end': gene_end,
            'feature': 'gene',
            'gene': gene_id
        })
        
        # 处理该基因的所有mRNA
        for mrna_info in mrnas[gene_id]:
            mrna_id = mrna_info['id']
            mrna_start = mrna_info['start']
            mrna_end = mrna_info['end']
            
            # 添加mRNA特征
            gene_dict[chrom].append({
                'start': mrna_start,
                'end': mrna_end,
                'feature': 'mRNA',
                'gene': gene_id
            })
            
            # 处理该mRNA的所有exon
            for exon in exons[mrna_id]:
                gene_dict[chrom].append({
                    'start': exon['start'],
                    'end': exon['end'],
                    'feature': 'exon',
                    'gene': gene_id
                })
            
            # 处理该mRNA的所有CDS
            for cds in cdss[mrna_id]:
                gene_dict[chrom].append({
                    'start': cds['start'],
                    'end': cds['end'],
                    'feature': 'CDS',
                    'gene': gene_id
                })
            
            # 计算intron区域 (mRNA中不属于exon的部分)
            if exons[mrna_id]:
                exons_sorted = sorted(exons[mrna_id], key=lambda x: x['start'])
                current_pos = mrna_start
                
                for exon in exons_sorted:
                    if exon['start'] > current_pos:
                        # 发现intron区域
                        gene_dict[chrom].append({
                            'start': current_pos,
                            'end': exon['start'],
                            'feature': 'intron',
                            'gene': gene_id
                        })
                    current_pos = max(current_pos, exon['end'])
                
                if current_pos < mrna_end:
                    # 最后一个intron
                    gene_dict[chrom].append({
                        'start': current_pos,
                        'end': mrna_end,
                        'feature': 'intron',
                        'gene': gene_id
                    })
    
    return gene_dict


def annotate_bedgraph(bedgraph_file, gff_dict):
    """注释bedGraph文件，按照特征优先级选择最佳匹配"""
    df = pd.read_csv(bedgraph_file, sep='\t', header=None, 
                     names=['chrom', 'start', 'end', 'signal'])
    
    # 定义特征优先级
    feature_priority = {
        'exon': 5,
        'intron': 4,
        'CDS': 3,
        'mRNA': 2,
        'gene': 1
    }
    
    gene_names = []
    features = []
    
    for _, row in df.iterrows():
        chrom = row['chrom']
        start = row['start']
        end = row['end']
        
        best_match = None
        best_priority = 0
        
        if chrom in gff_dict:
            for gene_info in gff_dict[chrom]:
                # 检查是否有重叠
                if not (end <= gene_info['start'] or start >= gene_info['end']):
                    feature = gene_info['feature']
                    priority = feature_priority.get(feature, 0)
                    
                    # 如果当前特征的优先级更高，则更新最佳匹配
                    if priority > best_priority:
                        best_match = gene_info
                        best_priority = priority
        
        if best_match:
            gene_names.append(best_match['gene'])
            features.append(best_match['feature'])
        else:
            gene_names.append('NA')
            features.append('Intergenic')
    
    df['gene_name'] = gene_names
    df['feature'] = features
    
    return df

def main():
    parser = argparse.ArgumentParser(description='Annotate bedGraph files with gene information')
    parser.add_argument('--bedgraph', required=True, help='Input bedGraph file')
    parser.add_argument('--gff', required=True, help='GFF annotation file')
    parser.add_argument('--output', required=True, help='Output CSV file')
    
    args = parser.parse_args()
    
    # 解析GFF文件
    gff_dict = parse_gff(args.gff)
    
    # 注释bedGraph文件
    annotated_df = annotate_bedgraph(args.bedgraph, gff_dict)
    
    # 保存结果
    annotated_df.to_csv(args.output, index=False)
    print(f"Annotated data saved to {args.output}")

if __name__ == '__main__':
    main()
