# 读取PLINK处理过的基因型数据（0/1/2），文件名为"data"
gdata <- fread("data.raw", header = T)

# 展示数据
gname <- colnames(gdata)
gdata[1:6, gname[1:10], with = F]

# 只取基因型
gdata <- gdata[ , !c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")]
gname <- colnames(gdata)
gdata[c(1:6), gname[1:6], with = F]