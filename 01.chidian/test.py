import pysam

# 查看pysam版本
print("pysam版本：", pysam.__version__)

# 验证BigWig接口是否存在
if hasattr(pysam, 'BigWig'):
    print("✅ 支持 pysam.BigWig（新版本）")
elif hasattr(pysam, 'BigWigFile'):
    print("✅ 支持 pysam.BigWigFile（旧版本）")
else:
    print("❌ 无BigWig相关接口，需重新安装pysam")
