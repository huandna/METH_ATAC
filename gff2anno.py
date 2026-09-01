import argparse
import pandas as pd
import re
def parse_gff_attributes(attributes_str):
    """解析GFF属性字段，提取ID和Name"""
    attributes = {}
    for item in attributes_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            attributes[key] = value
    return attributes
def convert_gff_to_anno(input_gff, output_anno):
    # 读取GFF文件
    gff_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gff_data = []
    with open(input_gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 9:
                gff_data.append(fields[:9])
    gff_df = pd.DataFrame(gff_data, columns=gff_columns)
    # 筛选基因特征（可根据需要调整类型）
    gene_df = gff_df[gff_df['type'] == 'gene'].copy()
    # 解析属性
    gene_df['attributes_dict'] = gene_df['attributes'].apply(parse_gff_attributes)
    gene_df['id'] = gene_df['attributes_dict'].apply(lambda x: x.get('ID', ''))
    gene_df['name'] = gene_df['attributes_dict'].apply(lambda x: x.get('Name', ''))
    # 构造输出格式
    anno_df = gene_df[['id', 'name', 'seqid', 'start', 'end']].copy()
    anno_df.columns = ['id', 'subanno', 'lineage', 'Sample', 'barcode']
    anno_df['anno'] = anno_df['id']  # 使用ID作为anno列
    anno_df = anno_df.fillna('')
    # 保存为制表符分隔文件
    anno_df.to_csv(output_anno, sep='\t', index=False, header=True)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert GFF file to annotation format.")
    parser.add_argument("--gff", required=True, help="Input GFF file path.")
    parser.add_argument("--output", required=True, help="Output annotation file path.")
    args = parser.parse_args()
    convert_gff_to_anno(args.gff, args.output)