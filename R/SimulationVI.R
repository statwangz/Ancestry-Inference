# 用变分推断方法处理模拟数据集

# 读取模拟数据集
g_data_s <- as.matrix(read.table("simulated_data.txt"))
theta_real <- as.matrix(read.table("theta.txt"))
N_s <- nrow(g_data_s)
L_s <- ncol(g_data_s)
K <- ncol(theta_real)

# 最大循环次数
MAX <- 10000

# 初始化
c <- 1/K
a <- 1
b <- 1
theta_s <- matrix(rgamma(N_s*K, 100, 0.01), nrow = N_s, ncol = K)
beta_s <- array(rbeta(K*L_s*2, a, b), dim = c(K, L_s, 2))
lower_bound <- numeric(length = MAX)

# 开始处理数据

# 将基因数据复制扩展为K*N*L维，便于后面的计算
g_data_temp <- array(t(apply(array(rep(g_data_s, K), dim = c(N_s, L_s, K)), MARGIN = 3, as.vector)), dim = c(K, N_s, L_s))

repeat_1 <- TRUE
s <- 1 # 计数器
lower_bound[s] <- LowerBound(g_data_s, theta_s, beta_s)

while(repeat_1){
  
  phi <- array(0, dim = c(K, N_s, L_s))
  xi <- array(0, dim = c(K, N_s, L_s))
  
  s <- s + 1
  
  beta_new <- array(dim = c(K, L_s, 2))
  phi_new <- array(dim = c(K, N_s, L_s))
  xi_new <- array(dim = c(K, N_s, L_s))
  
  # 对全数据更新 phi，xi 和 beta
  
  repeat_2 <- TRUE

  while(repeat_2){
    
    e_log_dir <- apply(theta_s, MARGIN = 1, FUN = ELogDir) # 输出结果每一列为每个人的 E_Log_Dir，K*N
    e_log_beta_1 <-  apply(beta_s, MARGIN = c(1, 2), FUN = ELogBeta1) # K*L
    e_log_beta_2 <-  apply(beta_s, MARGIN = c(1, 2), FUN = ELogBeta2)
    
    x <- array(as.vector(apply(e_log_beta_1, MARGIN = 2, FUN = function(x){x + e_log_dir})), dim = c(K, N_s, L_s))
    y <- array(as.vector(apply(e_log_beta_2, MARGIN = 2, FUN = function(x){x + e_log_dir})), dim = c(K, N_s, L_s))
    
    phi_new <- apply(x, MARGIN = c(2,3), FUN = LogSumExp) # K*N*L
    xi_new <- apply(y, MARGIN = c(2,3), FUN = LogSumExp)
    
    phi_temp <- g_data_temp * phi_new
    xi_temp <- (2 - g_data_temp) * xi_new
    
    beta_new[ , , 1] <- a + apply(phi_temp, MARGIN = c(1,3), sum) # K*L
    beta_new[ , , 2] <- b + apply(xi_temp, MARGIN = c(1,3), sum)
    
    if(Distance(beta_new, beta_s) + Distance(phi_new, phi) + Distance(xi_new, xi) < 1e-3){
      repeat_2 <- FALSE
    }
    
    phi <- phi_new
    xi <- xi_new
    beta_s <- beta_new
    
  }
  
  # 更新 theta
  theta_temp <- phi_temp + xi_temp
  theta_s <- c + t(apply(theta_temp, MARGIN = c(1,2), sum))
  
  # 计算 Lower Bound
  lower_bound[s] <- LowerBound(g_data_s, theta_s, beta_s)
  
  # 通过判断 Lower Bound 是否收敛来决定是否结束循环
  if((s == MAX) | (!is.na(lower_bound[s]) & !is.na(lower_bound[s-1]) & Distance(lower_bound[s], lower_bound[s-1]) < 1e-3)){
    repeat_1 <- FALSE
  }
  
}

E_theta <- matrix(nrow = N_s, ncol = K)
E_theta <- t(apply(theta_s, MARGIN = 1, FUN = function(x){x/sum(x)}))

# 结果与真实参数的差距
DMatrix(E_theta, theta_real)
write.table(E_theta, "theta_results_VI.txt", row.names = F, col.names = F)