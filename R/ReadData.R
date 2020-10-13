## 方法一：PLINK + data.table

# 读取PLINK处理过的基因型数据（0/1/2），文件名为"data"
g_data <- fread("data.raw", header = T)

# 展示数据
g_name <- colnames(g_data)
g_data[1:6, g_name[1:10], with = F]

# 只取基因型
g_data <- g_data[ , !c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")]
g_name <- colnames(g_data)
g_data[1:6, g_name[1:6], with = F]

## 方法二：genio

# 读取BED/BIM/FAM文件，文件名为"1kg_phase1_all_1m"
g_data_list <- read_plink("1kg_phase1_all_1m")
g_data <- as.data.table(t(g_data_list$X))
