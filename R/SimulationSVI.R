# 用随机变分推断方法处理模拟数据集

g_data_s <- as.matrix(read.table("simulated_data.txt"))
theta_s <- as.matrix(read.table("theta.txt"))

# 最大循环次数
MAX <- 10000

# 初始化
c <- 1/K
tau <- 1
kap <- 0.5
a <- 1
b <- 1
theta_s <- matrix(rgamma(N_s*K, 100, 0.01), nrow = N_s, ncol = K)
beta_s <- array(rbeta(K*L_s*2, a, b), dim = c(K, L_s, 2))
lower_bound <- numeric(length = MAX)

# 开始处理数据

repeat_1 <- TRUE
s <- 1 # 计数器
lower_bound[s] <- LowerBound(g_data_s, theta_s, beta_s)

while(repeat_1){
  
  phi <- matrix(0, nrow = N_s, ncol = K)
  xi <- matrix(0, nrow = N_s, ncol = K)
  
  s <- s + 1
  l <- sample(1:L_s, 1)
  g_data_l <- g_data_s[ , l]
  
  beta_new <- matrix(nrow = K, ncol = 2)
  phi_new <- matrix(nrow = N_s, ncol = K)
  xi_new <- matrix(nrow = N_s, ncol = K)
  
  repeat_2 <- TRUE
  
  while(repeat_2){
    
    e_log_dir <- apply(theta_s, MARGIN = 1, FUN = ELogDir)
    e_log_beta_1 <-  apply(beta_s[ , l, ], MARGIN = 1, FUN = ELogBeta1)
    e_log_beta_2 <-  apply(beta_s[ , l, ], MARGIN = 1, FUN = ELogBeta2)
    
    x <- e_log_dir + e_log_beta_1
    y <- e_log_dir + e_log_beta_2
    
    phi_new <- t(apply(x, MARGIN = 2, FUN = LogSumExp))
    xi_new <- t(apply(y, MARGIN = 2, FUN = LogSumExp))
    
    beta_new[ , 1] <- a + as.vector(t(g_data_l) %*% phi_new)
    beta_new[ , 2] <- b + as.vector(t(2 - g_data_l) %*% xi_new)
    
    if(Distance(beta_new, beta_s[ , l, ]) + Distance(phi_new, phi) + Distance(xi_new, xi) < 1e-3){
      repeat_2 <- FALSE
    }
    
    phi <- phi_new
    xi <- xi_new
    beta_s[ , l, ] <- beta_new
    
  }
  
  theta_new <- matrix(nrow = N_s, ncol = K)
  
  repeat_3 <- TRUE
  t <- 0
  
  while (repeat_3) {
    
    t <- t + 1
    
    rho <- (tau + t)^(-kap)
    theta_new <- (1 - rho) * theta_s + rho * (c + L_s * (g_data_l * phi + (2 - g_data_l) * xi))
    
    if(Distance(theta_new, theta_s) < 1e-3){
      repeat_3 <- FALSE
    }
    
    theta_s <- theta_new
    
  }
  
  # 计算 Lower Bound
  lower_bound[s] <- LowerBound(g_data_s, theta_s, beta_s)
  
  # 通过判断 Lower Bound 是否收敛来决定是否结束循环
  if((s == MAX) | (!is.na(lower_bound[s]) & !is.na(lower_bound[s-1]) & Distance(lower_bound[s], lower_bound[s-1]) < 1e-3)){
    repeat_1 <- FALSE
  }
  
}

# PCA
E_theta <- matrix(nrow = N_m, ncol = K)
E_theta <- t(apply(theta_s, MARGIN = 1, FUN = function(x){x/sum(x)}))

# 结果与真实参数的差距
Distance(E_theta, theta_real)
write.table(E_theta, "theta_results_SVI.txt", row.names = F, col.names = F)