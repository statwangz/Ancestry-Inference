repeat_1 <- TRUE
s <- 0 # 循环次数

theta_i <- array(dim = c(400, N, K)) # 存储 theta 迭代中的结果

while(repeat_1){
  
  s <- s + 1
  
  # 辅助参数
  phi <- matrix(0, nrow = N, ncol = K)
  xi <- matrix(0, nrow = N, ncol = K)
  
  # 抽取一列
  l <- sample(1:L, 1)
  g_data_l <- g_data[[l]]
  
  beta_new <- matrix(nrow = K, ncol = 2)
  phi_new <- matrix(nrow = N, ncol = K)
  xi_new <- matrix(nrow = N, ncol = K)
  
  repeat_2 <- TRUE
  
  while(repeat_2){
    
    e_log_dir <- apply(theta, MARGIN = 1, FUN = ELogDir) # 输出结果每一列为每个人的 E_Log_Dir
    e_log_beta_1 <-  apply(beta[ , l, ], MARGIN = 1, FUN = ELogBeta1)
    e_log_beta_2 <-  apply(beta[ , l, ], MARGIN = 1, FUN = ELogBeta2)
    
    x <- e_log_dir + e_log_beta_1 # 矩阵+向量，自动循环补齐，K*N，每一列为每个人的数据
    y <- e_log_dir + e_log_beta_2
    
    # log-sum-exp
    phi_new <- t(apply(x, MARGIN = 2, FUN = LogSumExp))
    xi_new <- t(apply(y, MARGIN = 2, FUN = LogSumExp))
    
    # 更新 beta
    beta_new[ , 1] <- a + as.vector(t(g_data_l) %*% phi_new)
    beta_new[ , 2] <- b + as.vector(t(2 - g_data_l) %*% xi_new)
    
    # 判断 phi, xi, beta 是否收敛
    if(Distance(beta_new, beta[ , l, ]) + Distance(phi_new, phi) + Distance(xi_new, xi) < 1e-3){
      repeat_2 <- FALSE
    }
    
    phi <- phi_new
    xi <- xi_new
    beta[ , l, ] <- beta_new
    
  }
  
  theta_new <- matrix(nrow = N, ncol = K)
  
  repeat_3 <- TRUE
  t <- 0 # 循环次数，用于计算 rho
  
  while (repeat_3) {
    
    t <- t + 1
    
    # 计算权重
    rho <- (tau + t)^(-kap)
    
    # 更新 theta
    theta_new <- (1 - rho) * theta + rho * (c + L * (g_data_l * phi + (2 - g_data_l) * xi))
    # 数*矩阵，向量*矩阵，自动循环补齐，N*K
    
    if(Distance(theta_new, theta) < 1e-3){
      repeat_3 <- FALSE
    }
    
    theta <- theta_new
    
  }
  
  # 每抽取 25 个基因存储一次 theta 值
  if((s%%25) == 0){
    theta_i[s/25, , ] <- theta
  }
  
  # 设定抽取基因个数，这里为 10000
  if(s == 10000){
    repeat_1 <- FALSE
  }
  
}