# 将分析真实数据得到的结果存储起来用于模拟
write.table(E_theta, "theta.txt", row.names = F, col.names = F)
write.table(E_beta[ , , 1], "beta.txt", row.names = F, col.names = F)

theta_m <- as.matrix(read.table("theta.txt"))
beta_m <- as.matrix(read.table("beta.txt"))

# 按照PSD model模拟基因数据

# 设置模拟数据集规模
N_m <- N    # 这里人数设置为与真实数据集相同
L_m <- 1000 # 这里基因数设置为1000

# 基因型0/1/2 ~ Binomial(2, p_il)，p_il = sum_k(theta_ik * beta_kl) 
p <- theta_m %*% beta_m
g_data_m <- matrix(rbinom(n = N_m*L_m, size = 2, prob = as.vector(p)), nrow = N_m, ncol = L_m)
write.table(g_data_m, "simulated_data.txt", row.names = F, col.names = F)
