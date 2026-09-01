#!/usr/bin/env python3
import subprocess
import argparse

def bedgraph_to_bigwig(bedgraph_file, bigwig_file, chrom_sizes_file):
    """将 bedGraph 转换为 BigWig"""
    cmd = [
        "bedGraphToBigWig",
        bedgraph_file,
        chrom_sizes_file,
        bigwig_file
    ]
    subprocess.run(cmd, check=True)

def main():
    parser = argparse.ArgumentParser(description="将 bedGraph 转换为 BigWig")
    parser.add_argument("-i", "--input", required=True, help="输入的 bedGraph 文件")
    parser.add_argument("-o", "--output", required=True, help="输出的 BigWig 文件")
    parser.add_argument("-c", "--chrom_sizes", required=True, help="染色体大小文件")
    args = parser.parse_args()

    bedgraph_to_bigwig(args.input, args.output, args.chrom_sizes)
    print(f"BigWig 文件已保存至: {args.output}")

if __name__ == "__main__":
    main()
