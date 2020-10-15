# 将分析真实数据得到的结果存储起来用于模拟
write.table(E_theta, "theta.txt", row.names = F, col.names = F)
write.table(E_beta[ , , 1], "beta.txt", row.names = F, col.names = F)

theta_s <- as.matrix(read.table("theta.txt"))
beta_s <- as.matrix(read.table("beta.txt"))

# 按照PSD model模拟基因数据

# 设置模拟数据集规模
N_s <- N    # 这里人数设置为与真实数据集相同
L_s <- 1000 # 这里基因数设置为1000

# 基因型0/1/2 ~ Binomial(2, p_il)，p_il = sum_k(theta_ik * beta_kl) 
p <- theta_s %*% beta_s
g_data_s <- as.data.table(matrix(rbinom(n = N_s*L_s, size = 2, prob = as.vector(p)), nrow = N_s, ncol = L_s))
fwrite(g_data_s, "simulated_data.txt", row.names = F, col.names = F)

# 下面开始验证

g_data_s <- fread("simulated_data.txt")

# 初始化

c <- 1/K
tau <- 1
kap <- 0.5
a <- 1
b <- 1
theta <- matrix(rgamma(N_s*K, 100, 0.01), nrow = N_s, ncol = K)
beta <- array(rbeta(K*L_s*2, a, b), dim = c(K, L_s, 2))

# 开始处理数据

repeat_1 <- TRUE
s <- 0 # 循环次数

while(repeat_1){
  
  # 辅助参数
  phi <- matrix(0, nrow = N_s, ncol = K)
  xi <- matrix(0, nrow = N_s, ncol = K)
  
  s <- s + 1
  
  g_data_l <- g_data_s[[s]]
  
  beta_new <- matrix(nrow = K, ncol = 2)
  phi_new <- matrix(nrow = N_s, ncol = K)
  xi_new <- matrix(nrow = N_s, ncol = K)
  
  repeat_2 <- TRUE
  
  while(repeat_2){
    
    e_log_dir <- apply(theta, MARGIN = 1, FUN = ELogDir) # 输出结果每一列为每个人的E_Log_Dir
    e_log_beta_1 <-  apply(beta[ , s, ], MARGIN = 1, FUN = ELogBeta1)
    e_log_beta_2 <-  apply(beta[ , s, ], MARGIN = 1, FUN = ELogBeta2)
    
    x <- e_log_dir + e_log_beta_1 # 矩阵+向量，自动循环补齐，K*N，每一列为每个人的数据
    y <- e_log_dir + e_log_beta_2
    
    # log-sum-exp
    phi_new <- t(apply(x, MARGIN = 2, FUN = LogSumExp))
    xi_new <- t(apply(y, MARGIN = 2, FUN = LogSumExp))
    
    # 更新beta参数
    beta_new[ , 1] <- a + as.vector(t(g_data_l) %*% phi_new)
    beta_new[ , 2] <- b + as.vector(t(2 - g_data_l) %*% xi_new)
    
    # 判断phi, xi, beta是否收敛
    if(DMatrix(beta_new, beta[ , s, ]) + DMatrix(phi_new, phi) + DMatrix(xi_new, xi) < 1e-3){
      repeat_2 <- FALSE
    }
    
    phi <- phi_new
    xi <- xi_new
    beta[ , s, ] <- beta_new
    
  }
 
  theta_new <- matrix(nrow = N_s, ncol = K)
  
  repeat_3 <- TRUE
  t <- 0 # 循环次数，用于计算rho
  
  while (repeat_3) {
    
    t <- t + 1
    # 计算权重
    rho <- (tau + t)^(-kap)
    
    # 更新theta
    theta_new <- (1 - rho) * theta + rho * (c + L_s * (g_data_l * phi + (2 - g_data_l) * xi))
    # 数*矩阵，向量*矩阵，自动循环补齐，N*K
    
    if(DMatrix(theta_new, theta) < 1e-3){
      repeat_3 <- FALSE
    }
    
    theta <- theta_new
    
  }
  
  if(s == L_m){
    repeat_1 <- FALSE
  }
  
}

# PCA
E_theta <- matrix(nrow = N_m, ncol = K)
E_theta <- t(apply(theta, MARGIN = 1, FUN = function(x){x/sum(x)}))
theta_PCA <- princomp(E_theta)
summary(theta_PCA, loadings = T)

PC_1 <- theta_PCA$scores[ , 1]
PC_2 <- theta_PCA$scores[ , 2]

pop_k4 <- pam(E_theta, k = 4)
table(pop_k4$clustering)

# 可视化
p <- ggplot(data = NULL, aes(x = PC_1, y = PC_2))
p + geom_point(color = pop_k4$clustering + 1, alpha = 1/2, size = 2) +
  scale_x_continuous(name = "PC 1", limits = range(PC_1)) +
  scale_y_continuous(name = "PC 2", limits = range(PC_2)) +
  labs(title="PCA")

# 结果与真实参数的差距
DMatrix(E_theta, theta_s)