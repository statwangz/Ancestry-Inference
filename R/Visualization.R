# 查看哪次迭代后的结果，这里选25*12
i <- 12

# 计算theta期望
E_theta <- matrix(nrow = N, ncol = K)
E_theta <- t(apply(theta_r[i, , ], MARGIN = 1, FUN = function(x){x/sum(x)}))

# PCA
theta_PCA <- princomp(E_theta)
summary(theta_PCA, loadings = T)
PC_1 <- theta_PCA$scores[ , 1]
PC_2 <- theta_PCA$scores[ , 2]

# 聚类，这里选定聚类数为4
pop_k4 <- pam(E_theta, k = 4)
table(pop_k4$clustering)

# 绘图
ggplot(data = NULL, aes(x = PC_1, y = PC_2)) +
  geom_point(color = pop_k4$cluster + 1) +
  labs(title="PCA", x = "PC 1", y = "PC 2")