# 用于判断是否收敛
DMatrix <- function(matrix_new, matrix_old){
  d <- as.vector(matrix_new - matrix_old)
  sum(abs(d))/sum(abs(matrix_new))
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