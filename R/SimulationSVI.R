# 用随机变分推断方法处理模拟数据集

g_data_s <- as.matrix(read.table("simulated_data.txt"))
theta_real <- as.matrix(read.table("theta_real.txt"))
N_s <- nrow(g_data_s)
L_s <- ncol(g_data_s)
K <- ncol(theta_real)

# 最大循环次数
MAX <- floor(L_s/4) # 取总基因数的 25% 进行计算

# 初始化
c <- 1/K
tau <- 1
kap <- 0.5
a <- 1
b <- 1
theta_s <- matrix(rgamma(N_s*K, 100, 0.01), nrow = N_s, ncol = K)
beta_s <- array(rbeta(K*L_s*2, a, b), dim = c(K, L_s, 2))

# 储存 Lower Bound，每更新一列基因数据更新一次
# 一般情况下当抽取基因比例较小时，Lower Bound 很小（考虑没有参与计算的 beta），所以没有计算的意义
lower_bound <- numeric(length = MAX)

# 储存更新后的 theta， 每更新 50 次储存一次
theta_array <- array(dim = c(floor(MAX/50), N_s, K))

# 开始处理数据

repeat_1 <- TRUE # 控制循环
s <- 0 # 计数器
lower_bound_0 <- LowerBound(g_data_s, theta_s, beta_s)

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
  
  lower_bound[s] <- LowerBound(g_data_s, theta_s, beta_s)

  # 记录更新过程中的 theta
  if(s%%50 == 0){
    theta_array[s/50, , ] <- theta_s
  }
  
  # 更新 MAX 次后停止
  if(s == MAX){
    repeat_1 <- FALSE
  }
  
}

# 计算期望
E_theta <- matrix(nrow = N_s, ncol = K)
E_theta <- t(apply(theta_s, MARGIN = 1, FUN = function(x){x/sum(x)}))

# 结果与真实参数的差距
Distance(E_theta, theta_real)
write.table(E_theta, "theta_results_SVI.txt", row.names = F, col.names = F)