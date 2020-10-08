repeat_1 <- TRUE
s <- 0 # 循环次数

theta_i <- array(dim = c(400, N, K)) # 存储theta迭代中的结果

while(repeat_1){
  
  s <- s + 1
  
  # 辅助参数
  phi <- matrix(0, nrow = N, ncol = K)
  xi <- matrix(0, nrow = N, ncol = K)
  
  l <- sample(1:L, 1)
  gdata_l <- gdata[[l]]
  
  beta_new <- matrix(nrow = K, ncol = 2)
  phi_new <- matrix(nrow = N, ncol = K)
  xi_new <- matrix(nrow = N, ncol = K)
  
  repeat_2 <- TRUE
  
  while(repeat_2){
    
    e_log_dir <- apply(theta, MARGIN = 1, FUN = ELogDir) # 输出结果每一列为每个人的E_Log_Dir
    e_log_beta_1 <-  apply(beta[ , l, ], MARGIN = 1, FUN = ELogBeta1)
    e_log_beta_2 <-  apply(beta[ , l, ], MARGIN = 1, FUN = ELogBeta2)
    
    x <- e_log_dir + e_log_beta_1 # 矩阵+向量，自动循环补齐，K*N，每一列为每个人的数据
    y <- e_log_dir + e_log_beta_2
    
    # log-sum-exp
    phi_new <- t(apply(x, MARGIN = 2, FUN = LogSumExp))
    xi_new <- t(apply(y, MARGIN = 2, FUN = LogSumExp))
    
    # 更新beta参数
    beta_new[ , 1] <- a + as.vector(t(gdata_l) %*% phi_new)
    beta_new[ , 2] <- b + as.vector(t(2 - gdata_l) %*% xi_new)
    
    # 判断phi, xi, beta是否收敛
    if(DMatrix(beta_new, beta[ , l, ]) + DMatrix(phi_new, phi) + DMatrix(xi_new, xi) < 1e-3){
      repeat_2 <- FALSE
    }
    
    phi <- phi_new
    xi <- xi_new
    beta[ , l, ] <- beta_new
    
  }
  
  theta_new <- matrix(nrow = N, ncol = K)
  
  repeat_3 <- TRUE
  t <- 0 # 循环次数，用于计算rho
  
  while (repeat_3) {
    
    t <- t + 1
    
    # 计算权重
    rho <- (tau + t)^(-kap)
    
    # 更新theta
    theta_new <- (1 - rho) * theta + rho * (c + L * (gdata_l * phi + (2 - gdata_l) * xi))
    # 数*矩阵，向量*矩阵，自动循环补齐，N*K
    
    if(DMatrix(theta_new, theta) < 1e-3){
      repeat_3 <- FALSE
    }
    
    theta <- theta_new
    
  }
  
  # 每抽取25个基因存储一次theta值
  if((s%%25) == 0){
    theta_i[s/25, , ] <- theta
  }
  
  # 设定抽取基因个数，这里为10000
  if(s == 10000){
    repeat_1 <- FALSE
  }
  
}