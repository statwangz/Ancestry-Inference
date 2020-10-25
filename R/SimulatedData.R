# 将分析真实数据得到的结果存储起来用于模拟
write.table(E_theta, "theta.txt", row.names = F, col.names = F)
write.table(E_beta[ , , 1], "beta.txt", row.names = F, col.names = F)

theta_real <- as.matrix(read.table("theta.txt")) # 真实血统比例
beta_real <- as.matrix(read.table("beta.txt"))   # 真实基因频率

# 按照 PSD model 模拟基因数据

# 设置模拟数据集规模
N_s <- N    # 这里人数设置为与真实数据集相同
L_s <- 1000 # 这里基因数设置为 1000
K <- 3      # K 依然为3 

# 基因型 X = 0/1/2 ~ Binomial(2, p_il)，p_il = sum_k(theta_ik * beta_kl)
p <- theta_real %*% beta_real
g_data_s <- as.data.table(matrix(rbinom(n = N_s*L_s, size = 2, prob = as.vector(p)), nrow = N_s, ncol = L_s))
fwrite(g_data_s, "simulated_data.txt", row.names = F, col.names = F)

# 下面开始验证

g_data_s <- as.matrix(read.table("simulated_data.txt"))

# 初始化

c <- 1/K
tau <- 1
kap <- 0.5
a <- 1
b <- 1
theta_s <- matrix(rgamma(N_s*K, 100, 0.01), nrow = N_s, ncol = K)
beta_s <- array(rbeta(K*L_s*2, a, b), dim = c(K, L_s, 2))

# 开始处理数据

repeat_1 <- TRUE
s <- 0

lower_bound <- numeric(length = L_s)

while(repeat_1){
  
  phi <- matrix(0, nrow = N_s, ncol = K)
  xi <- matrix(0, nrow = N_s, ncol = K)
  
  s <- s + 1
  
  g_data_l <- g_data_s[[s]]
  
  beta_new <- matrix(nrow = K, ncol = 2)
  phi_new <- matrix(nrow = N_s, ncol = K)
  xi_new <- matrix(nrow = N_s, ncol = K)
  
  repeat_2 <- TRUE
  
  while(repeat_2){
    
    e_log_dir <- apply(theta_s, MARGIN = 1, FUN = ELogDir)
    e_log_beta_1 <-  apply(beta_s[ , s, ], MARGIN = 1, FUN = ELogBeta1)
    e_log_beta_2 <-  apply(beta_s[ , s, ], MARGIN = 1, FUN = ELogBeta2)
    
    x <- e_log_dir + e_log_beta_1
    y <- e_log_dir + e_log_beta_2
    
    phi_new <- t(apply(x, MARGIN = 2, FUN = LogSumExp))
    xi_new <- t(apply(y, MARGIN = 2, FUN = LogSumExp))
    
    beta_new[ , 1] <- a + as.vector(t(g_data_l) %*% phi_new)
    beta_new[ , 2] <- b + as.vector(t(2 - g_data_l) %*% xi_new)
    
    if(DMatrix(beta_new, beta_s[ , s, ]) + DMatrix(phi_new, phi) + DMatrix(xi_new, xi) < 1e-3){
      repeat_2 <- FALSE
    }
    
    phi <- phi_new
    xi <- xi_new
    beta_s[ , s, ] <- beta_new
    
  }
 
  theta_new <- matrix(nrow = N_s, ncol = K)
  
  repeat_3 <- TRUE
  t <- 0
  
  while (repeat_3) {
    
    t <- t + 1
    
    rho <- (tau + t)^(-kap)
    theta_new <- (1 - rho) * theta_s + rho * (c + L_s * (g_data_l * phi + (2 - g_data_l) * xi))
    
    if(DMatrix(theta_new, theta_s) < 1e-3){
      repeat_3 <- FALSE
    }
    
    theta_s <- theta_new
    
  }
  
  if(s == L_s){
    repeat_1 <- FALSE
  }
  
  lower_bound[s] <- LowerBound(g_data_s, theta_s, beta_s, K, N_s, L_s)
  
}

# PCA
E_theta <- matrix(nrow = N_m, ncol = K)
E_theta <- t(apply(theta_s, MARGIN = 1, FUN = function(x){x/sum(x)}))
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
DMatrix(E_theta, theta_real)
write.table(E_theta, "theta_results.txt", row.names = F, col.names = F)