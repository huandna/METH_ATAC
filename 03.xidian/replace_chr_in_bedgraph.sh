#!/bin/bash
# ============================================================
# bedGraph 染色体编号替换脚本 (Linux)
# 读取 chrs_paired.xls (TSV: 旧编号<TAB>新编号)
# 替换 bedGraph 第1列旧编号→新编号，不在映射表中的行舍弃
#
# 用法:
#   bash replace_chr_in_bedgraph.sh <input.bedGraph> [output.bedGraph] [chrs_paired.xls]
#
# 示例:
#   bash replace_chr_in_bedgraph.sh sample.bedGraph
#   bash replace_chr_in_bedgraph.sh sample.bedGraph output_renamed.bedGraph
#   bash replace_chr_in_bedgraph.sh sample.bedGraph out.bg /path/to/chrs_paired.xls
# ============================================================

set -euo pipefail

# ---- 参数解析 ----
INPUT="${1:?用法: bash replace_chr_in_bedgraph.sh <input.bedGraph> [output.bedGraph] [chrs_paired.xls]}"
OUTPUT="${2:-${INPUT%.*}_chrRenamed.bedGraph}"

# 映射表路径: 优先用第3参数 → 同级目录 → 项目根目录
CHRS_MAP="${3:-}"
if [ -z "$CHRS_MAP" ]; then
    SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
    if [ -f "${SCRIPT_DIR}/chrs_paired.xls" ]; then
        CHRS_MAP="${SCRIPT_DIR}/chrs_paired.xls"
    elif [ -f "${SCRIPT_DIR}/../chrs_paired.xls" ]; then
        CHRS_MAP="${SCRIPT_DIR}/../chrs_paired.xls"
    else
        CHRS_MAP="chrs_paired.xls"
    fi
fi

# ---- 检查文件 ----
if [ ! -f "$INPUT" ]; then
    echo "[错误] 输入文件不存在: $INPUT" >&2; exit 1
fi
if [ ! -f "$CHRS_MAP" ]; then
    echo "[错误] 映射文件不存在: $CHRS_MAP" >&2; exit 1
fi

echo "=========================================="
echo "bedGraph 染色体编号替换"
echo "=========================================="
echo "  输入: $INPUT"
echo "  输出: $OUTPUT"
echo "  映射: $CHRS_MAP"
echo ""

# ---- awk 一趟完成: 读映射 + 替换/过滤 ----
awk -F'\t' -v OFS='\t' '
# 第一个文件: 映射表 → 关联数组
NR == FNR {
    chr_map[$1] = $2
    next
}
# 第二个文件: bedGraph
{
    # 保留 header 行 (browser / track 开头)
    if ($1 == "browser" || $1 == "track") {
        print $0
        next
    }
    # 第1列查找映射
    if ($1 in chr_map) {
        $1 = chr_map[$1]
        print $0
        kept++
    } else {
        skipped++
    }
}
END {
    printf "[完成] 保留: %d 行, 舍弃: %d 行\n", kept+0, skipped+0 > "/dev/stderr"
}
' "$CHRS_MAP" "$INPUT" > "$OUTPUT"

echo ""
echo "输出: $OUTPUT ($(wc -l < "$OUTPUT") 行)"
