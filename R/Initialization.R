K <- 3            # 预设祖先群体数，这里设定为 3
N <- nrow(g_data) # 人数
L <- ncol(g_data) # 检测位点数

## 超参数设定

# Dirichlet
c <- 1/K

# 用于计算 theta 参数更新时的权重
tau <- 1
kap <- 0.5

# 初始 Beta 参数
a <- 1
b <- 1

# 初始化 theta (每个人的血统比例的分布参数，N*K，Dirichlet-多项式分布)
theta <- matrix(rgamma(N*K, 100, 0.01), nrow = N, ncol = K) # 从 gamma(100, 0.01) 抽取

# 初始化 beta (每个群体在每个位点上的基因型的分布参数，K*L*2，Beta-二项分布)
beta <- array(rbeta(K*L*2, a, b), dim = c(K, L, 2)) # 从 beta(a,b) 抽取