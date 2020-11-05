# 查看哪次迭代后的结果，这里选 25*12
i <- 12

# 计算 theta 期望
E_theta <- matrix(nrow = N, ncol = K)
E_theta <- t(apply(theta_r[i, , ], MARGIN = 1, FUN = function(x){x/sum(x)}))

# PCA
theta_PCA <- princomp(E_theta)
summary(theta_PCA, loadings = T)
PC_1 <- theta_PCA$scores[ , 1]
PC_2 <- theta_PCA$scores[ , 2]

# 聚类，这里选定聚类数为 4
pop_k4 <- pam(E_theta, k = 4)
table(pop_k4$clustering)

# 绘图
p <- ggplot(data = NULL, aes(x = PC_1, y = PC_2))
p + geom_point(color = pop_k4$clustering + 1, alpha = 1/2, size = 2) +
  scale_x_continuous(name = "PC 1", limits = range(PC_1)) +
  scale_y_continuous(name = "PC 2", limits = range(PC_2)) +
  labs(title="PCA")

# 血统比例图
theta_plot <- cbind(c(1:N), E_theta)
colnames(theta_plot) <- c("ID", "K1", "K2", "K3")
theta_plot <- melt(as.data.frame(theta_plot), id.vars = "ID", variable.name = "K", value.name = "Proportion")
ggplot(theta_plot, aes(ID, Proportion, fill = K)) +
  geom_bar(stat = "identity", position = "fill")