# 用于判断是否收敛
Distance <- function(new, old){
  d <- abs(new - old)
  sum(d)/sum(abs(new))
}

# X ~ Beta，计算 ln(X) 的期望
ELogBeta1 <- function(beta){
  digamma(beta[1]) - digamma(sum(beta))
}
# X ~ Beta，计算 ln(1-X) 的期望
ELogBeta2 <- function(beta){
  digamma(beta[2]) - digamma(sum(beta))
}
# X ~ Dirichlet，计算 ln(X_k) 的期望
ELogDir <- function(dir){
  digamma(dir) - digamma(sum(dir))
}

# log-sum-exp
LogSumExp <- function(x){
  x <- x - max(x)
  exp(x - log(sum(exp(x))))
}

# 计算 Lower Bound
LowerBound <- function(g_data, theta, beta){
  
  N <- nrow(g_data)
  L <- ncol(g_data)
  K <- ncol(theta)
  
  e_log_dir <- apply(theta, MARGIN = 1, FUN = ELogDir) # 输出结果每一列为每个人的 E_Log_Dir，K*N
  e_log_beta_1 <-  apply(beta, MARGIN = c(1, 2), FUN = ELogBeta1) # K*L
  e_log_beta_2 <-  apply(beta, MARGIN = c(1, 2), FUN = ELogBeta2)
  
  E_1 <- sum(as.vector(a - beta[ , , 1]) * as.vector(e_log_beta_1))
  
  E_2 <- sum(as.vector(b - beta[ , , 2]) * as.vector(e_log_beta_2))
  
  E_3 <- sum(as.vector(t(c - theta)) * as.vector(e_log_dir))
  
  # 输出结果每一列为每个基因的信息，相当于对每个基因抽样分析后得到的x矩阵
  x <- array(as.vector(apply(e_log_beta_1, MARGIN = 2, FUN = function(x){x + e_log_dir})), dim = c(K, N, L))
  y <- array(as.vector(apply(e_log_beta_2, MARGIN = 2, FUN = function(x){x + e_log_dir})), dim = c(K, N, L))
  
  phi <- apply(x, MARGIN = c(2,3), FUN = LogSumExp) # K*N*L
  xi <- apply(y, MARGIN = c(2,3), FUN = LogSumExp)
  
  e_log_dir_temp <- array(rep(e_log_dir, L), dim = c(K, N, L))
  e_log_beta_1_temp <- array(t(apply(array(rep(e_log_beta_1, N), dim = c(K, L, N)), MARGIN = 1, t)), dim = c(K, N, L))
  e_log_beta_2_temp <- array(t(apply(array(rep(e_log_beta_2, N), dim = c(K, L, N)), MARGIN = 1, t)), dim = c(K, N, L))
  
  
  E_4 <- sum(as.vector(g_data) * as.vector(apply(array(as.vector(phi) * as.vector(e_log_dir_temp + e_log_beta_1_temp - log(phi)), dim = c(K, N, L)), MARGIN = c(2, 3), sum)))
  E_5 <- sum(as.vector(2 - g_data) * as.vector(apply(array(as.vector(xi) * as.vector(e_log_dir_temp + e_log_beta_2_temp - log(xi)), dim = c(K, N, L)), MARGIN = c(2, 3), sum)))
  
  return(E_1 + E_2 + E_3 + E_4 + E_5)
  
}