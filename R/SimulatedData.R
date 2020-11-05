### 采用两种方法生成模拟数据集

## 数据集一：将分析真实数据得到的结果用于模拟

# 设置数据集规模：人数与真实数据集相同，基因选择用 SVI 方法抽取的前 1000 个
# N_s = N，人数设置为与真实数据集相同
# L_s = 1000，基因数设置为 1000
# K = 3，K 依然为3 

# 存储
write.table(E_theta, "theta_real.txt", row.names = F, col.names = F) # N*K
write.table(E_beta, "beta_real.txt", row.names = F, col.names = F)   # K*L

# 读取
theta_real <- as.matrix(read.table("theta_real.txt")) # 真实血统比例
beta_real <- as.matrix(read.table("beta_real.txt"))   # 真实基因频率

N_s <- nrow(theta_real)
L_s <- ncol(beta_real)
K <- nrow(beta_real)

# 基因型 X = 0/1/2 ~ Binomial(2, p_il)，p_il = sum_k(theta_ik * beta_kl)
p <- theta_real %*% beta_real
g_data_s <- as.data.table(matrix(rbinom(n = N_s*L_s, size = 2, prob = as.vector(p)), nrow = N_s, ncol = L_s))
write.table(g_data_s, "simulated_data.txt", row.names = F, col.names = F)

## 数据集二：即原 paper 中生成模拟数据集的第二种方法

# 基因频率与一中相同，血统比例采用混合的高斯分布得到，具有空间相关性
# 设置数据集规模：人数为1000，基因选择用 SVI 方法抽取的前 2000 个
# N_s = 1000
# L_s = 2000
# K = 3，K 依然为3 

# 存储基因频率
write.table(E_beta, "beta_real.txt", row.names = F, col.names = F) # K*L

# 读取基因频率
beta_real <- as.matrix(read.table("beta_real.txt")) # 真实基因频率
L_s <- ncol(beta_real)
K <- nrow(beta_real)

# 生成血统比例
N_s <- 1000
x <- (K + 1)/(N_s - 1) * c(0 : (N_s-1)) # 在区间 [0, K+1] 上等间距生成 N_s 个点，即为每个人坐标
mu <- c(1:K)   # 设置高斯分布均值，即整数点坐标
s <- rep(1, K) # 设置高斯分布方差，这里都取 2
theta_temp <- matrix(mapply(dnorm, rep(x, K), rep(mu, each = N_s), rep(s, each = N_s)), nrow = N_s, ncol = K)
theta_real <- t(apply(theta_temp, MARGIN = 1, FUN = function(x){x/sum(x)})) # 归一化

# 存储血统比例
write.table(theta_real, "theta_real.txt", row.names = F, col.names = F) # N*K

# 生成基因型并存储，同方法一
p <- theta_real %*% beta_real
g_data_s <- as.data.table(matrix(rbinom(n = N_s*L_s, size = 2, prob = as.vector(p)), nrow = N_s, ncol = L_s))
write.table(g_data_s, "simulated_data.txt", row.names = F, col.names = F)