# 将分析真实数据得到的结果存储起来用于模拟
write.table(E_theta, "theta.txt", row.names = F, col.names = F)
write.table(E_beta[ , , 1], "beta.txt", row.names = F, col.names = F)

theta_real <- as.matrix(read.table("theta.txt")) # 真实血统比例
beta_real <- as.matrix(read.table("beta.txt"))   # 真实基因频率

## 按照 PSD model 模拟基因数据
# 设置模拟数据集规模
N_s <- N    # 这里人数设置为与真实数据集相同
L_s <- 1000 # 这里基因数设置为 1000
K <- 3      # K 依然为3 
# 基因型 X = 0/1/2 ~ Binomial(2, p_il)，p_il = sum_k(theta_ik * beta_kl)
p <- theta_real %*% beta_real
g_data_s <- as.data.table(matrix(rbinom(n = N_s*L_s, size = 2, prob = as.vector(p)), nrow = N_s, ncol = L_s))
fwrite(g_data_s, "simulated_data.txt", row.names = F, col.names = F)