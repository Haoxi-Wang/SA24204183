## -----------------------------------------------------------------------------
library(ggplot2)

# 自定义函数：生成雷利分布的随机样本
# n: 样本量
# sigma: 分布的参数 σ（理论众数）
generate_rayleigh <- function(n, sigma) {
  # 使用逆变换法生成雷利分布的随机样本
  # U 是 (0, 1) 上的均匀分布
  U <- runif(n)
  # 使用雷利分布的逆CDF公式生成样本
  rayleigh_samples <- sigma * sqrt(-2 * log(1 - U))
  
  return(rayleigh_samples)
}

# 自定义函数：绘制直方图并检查众数
plot_rayleigh_histogram <- function(samples, sigma) {
  # 将生成的样本转换为数据框，以便使用 ggplot 绘图
  data <- data.frame(samples = samples)
  
  # 绘制直方图
  p <- ggplot(data, aes(x = samples)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, color = "black", fill = "skyblue", alpha = 0.7) +
    # 添加理论雷利分布的密度曲线
    stat_function(fun = function(x) x / sigma^2 * exp(-x^2 / (2 * sigma^2)), color = "red", linewidth = 1) +
    # 添加理论众数 σ 的垂直线
    geom_vline(xintercept = sigma, color = "blue", linetype = "dashed", linewidth = 1) +
    labs(title = paste("Histogram of Rayleigh Samples (sigma =", sigma, ")"),
         x = "Samples", y = "Density") +
    theme_minimal()
  
  print(p)
}

# 设置样本量和 sigma 的不同值
n <- 10000  # 样本量
sigma_values <- c(1, 2, 3)  # 不同的 σ 值

# 对每个 σ 值生成样本并绘制直方图
for (sigma in sigma_values) {
  samples <- generate_rayleigh(n, sigma)  # 生成样本
  plot_rayleigh_histogram(samples, sigma)  # 绘制直方图并检查众数
}


## -----------------------------------------------------------------------------
# 设置随机种子，以便结果可重复
set.seed(123)

# 混合分布的参数
n <- 1000        # 样本大小
p1 <- 0.75       # N(0, 1) 的混合概率
p2 <- 1 - p1     # N(3, 1) 的混合概率

# 生成混合正态分布的随机样本
# 先生成一个大小为n的样本，确定每个样本是来自哪一个分布
# rbinom(n, 1, p1) 返回一个0和1的向量，1代表来自N(0,1)，0代表来自N(3,1)
mix_sample <- rbinom(n, 1, p1)

# 根据 mix_sample 的值生成最终的样本
# 如果是1，则从N(0, 1)生成值；如果是0，则从N(3, 1)生成值
sample <- ifelse(mix_sample == 1, rnorm(n, 0, 1), rnorm(n, 3, 1))

# 创建数据框，方便使用 ggplot
df <- data.frame(sample = sample)

# 绘制直方图，叠加核密度估计曲线
ggplot(df, aes(x = sample)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, color = "black", fill = "lightblue") +  # 直方图
  geom_density(color = "red", linewidth = 1.2) +  # 核密度曲线
  ggtitle("Normal Mixture Distribution with p1 = 0.75") +  # 图标题
  xlab("Value") +  # X轴标签
  ylab("Density") +  # Y轴标签
  theme_minimal()  # 简洁的主题

# 为不同的 p1 值绘制分布，观察双峰现象
# p1 值可以设置为 0.5、0.25、0.1 等
p_values <- c(0.5, 0.25, 0.1)

for (p1 in p_values) {
  p2 <- 1 - p1
  mix_sample <- rbinom(n, 1, p1)
  sample <- ifelse(mix_sample == 1, rnorm(n, 0, 1), rnorm(n, 3, 1))
  
  df <- data.frame(sample = sample)
  
  # 绘制每个p1值的直方图和密度曲线
  p <- ggplot(df, aes(x = sample)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, color = "black", fill = "lightblue") +
    geom_density(color = "red", linewidth = 1.2) +
    ggtitle(paste("Normal Mixture Distribution with p1 =", p1)) +
    xlab("Value") +
    ylab("Density") +
    theme_minimal()
  # 打印绘图
  print(p)
}

## -----------------------------------------------------------------------------
# 设置随机种子以确保可重复性
set.seed(123)

# 泊松-伽马过程参数
lambda <- 2  # 泊松过程参数 λ
t <- 10      # 时间 t
alpha <- 2   # 伽马分布形状参数
beta <- 1    # 伽马分布尺度参数

# 模拟复合泊松过程
simulate_compound_poisson_gamma <- function(lambda, t, alpha, beta) {
  # 生成泊松过程N(t)
  N_t <- rpois(1, lambda * t)
  
  # 如果N(t) > 0，生成N(t)个伽马随机变量Y_i
  if (N_t > 0) {
    Y <- rgamma(N_t, shape = alpha, scale = beta)
    X_t <- sum(Y)  # 复合泊松过程 X(t)
  } else {
    X_t <- 0
  }
  
  return(X_t)
}

# 模拟多个复合泊松过程并估计均值和方差
n_simulations <- 10000  # 模拟次数
X_t_samples <- replicate(n_simulations, simulate_compound_poisson_gamma(lambda, t, alpha, beta))

# 估计的均值和方差
estimated_mean <- mean(X_t_samples)
estimated_variance <- var(X_t_samples)

# 理论值计算
theoretical_mean <- lambda * t * alpha * beta
theoretical_variance <- lambda * t * (alpha * beta^2 + alpha^2 * beta^2)

# 输出结果
cat("估计均值: ", estimated_mean, "\n")
cat("理论均值: ", theoretical_mean, "\n")
cat("估计方差: ", estimated_variance, "\n")
cat("理论方差: ", theoretical_variance, "\n")

#更换不同的多组参数值，比较模拟结果与理论值
lambda <- 3  # 泊松过程参数 λ
t <- 10      # 时间 t
alpha <- 3   # 伽马分布形状参数
beta <- 2    # 伽马分布尺度参数
X_t_samples <- replicate(n_simulations, simulate_compound_poisson_gamma(lambda, t, alpha, beta))
estimated_mean <- mean(X_t_samples)
estimated_variance <- var(X_t_samples)
theoretical_mean <- lambda * t * alpha * beta
theoretical_variance <- lambda * t * (alpha * beta^2 + alpha^2 * beta^2)
cat("估计均值: ", estimated_mean, "\n")
cat("理论均值: ", theoretical_mean, "\n")
cat("估计方差: ", estimated_variance, "\n")
cat("理论方差: ", theoretical_variance, "\n")

lambda <- 4  # 泊松过程参数 λ
t <- 10      # 时间 t
alpha <- 4   # 伽马分布形状参数
beta <- 3    # 伽马分布尺度参数
X_t_samples <- replicate(n_simulations, simulate_compound_poisson_gamma(lambda, t, alpha, beta))
estimated_mean <- mean(X_t_samples)
estimated_variance <- var(X_t_samples)
theoretical_mean <- lambda * t * alpha * beta
theoretical_variance <- lambda * t * (alpha * beta^2 + alpha^2 * beta^2)
cat("估计均值: ", estimated_mean, "\n")
cat("理论均值: ", theoretical_mean, "\n")
cat("估计方差: ", estimated_variance, "\n")
cat("理论方差: ", theoretical_variance, "\n")

## -----------------------------------------------------------------------------
# 设置随机种子，确保结果可重复
set.seed(123)

# 定义函数 MonteCarloBeta，用于计算 Beta(3,3) 分布的 CDF 的蒙特卡洛估计
MonteCarloBeta <- function(x, n_sim = 10000) {
  # 从 Beta(3,3) 分布中生成 n_sim 个随机样本
  beta_samples <- rbeta(n_sim, shape1 = 3, shape2 = 3)
  
  # 计算 F(x) 的蒙特卡洛估计，即随机样本中小于等于 x 的比例
  F_x_hat <- mean(beta_samples <= x)
  
  return(F_x_hat)
}

# 要计算的 x 值
x_values <- seq(0.1, 0.9, by = 0.1)

# 初始化一个空的数据框来存储结果
results <- data.frame(x = x_values, MonteCarlo = NA, pbeta = NA, Difference = NA)

# 进行 Monte Carlo 估计和使用 pbeta 函数的比较
for (i in 1:length(x_values)) {
  # 计算 Monte Carlo 估计
  results$MonteCarlo[i] <- MonteCarloBeta(x_values[i])
  
  # 使用 pbeta 函数计算 Beta(3,3) 的精确 CDF 值
  results$pbeta[i] <- pbeta(x_values[i], shape1 = 3, shape2 = 3)
  
  # 计算 Monte Carlo 估计与 pbeta 值的差异
  results$Difference[i] <- abs(results$MonteCarlo[i] - results$pbeta[i])
}

# 打印结果
print(results)

# 画图比较 Monte Carlo 估计和 pbeta 函数的值
plot(x_values, results$MonteCarlo, type = "b", col = "blue", pch = 16, ylim = c(0, 1),
     xlab = "x", ylab = "CDF", main = "Monte Carlo vs. pbeta CDF Estimation")
lines(x_values, results$pbeta, type = "b", col = "red", pch = 16)
legend("bottomright", legend = c("Monte Carlo Estimate", "pbeta Value"), col = c("blue", "red"), lty = 1, pch = 16)

## -----------------------------------------------------------------------------
# 设置随机种子以确保可重复性
set.seed(123)

# Rayleigh分布的随机变量生成函数
# 生成Rayleigh(σ)分布的样本
rayleigh_sample <- function(n, sigma = 1) {
  # 标准正态分布样本
  z <- rnorm(n, mean = 0, sd = sigma)
  
  # Rayleigh分布样本（z^2 开方）
  x <- sqrt(z^2)
  
  return(x)
}

# 生成使用对偶变量的Rayleigh样本
# 使用对偶变量生成样本并计算平均值
rayleigh_antithetic <- function(n, sigma = 1) {
  # 生成均匀分布样本
  u <- runif(n)
  
  # 使用对偶变量：生成 u 和 1-u
  x1 <- sqrt(-2 * sigma^2 * log(u))
  x2 <- sqrt(-2 * sigma^2 * log(1 - u))
  
  # 对偶变量的均值
  avg_antithetic <- (x1 + x2) / 2
  
  return(avg_antithetic)
}

# 生成独立Rayleigh样本并计算均值
rayleigh_independent <- function(n, sigma = 1) {
  x1 <- sqrt(-2 * sigma^2 * log(runif(n)))
  x2 <- sqrt(-2 * sigma^2 * log(runif(n)))
  
  # 独立样本的均值
  avg_independent <- (x1 + x2) / 2
  
  return(avg_independent)
}

# 样本数量
n <- 10000
sigma <- 1

# 计算使用对偶变量的样本均值
antithetic_samples <- rayleigh_antithetic(n, sigma)
# 计算独立样本的均值
independent_samples <- rayleigh_independent(n, sigma)

# 计算方差
var_antithetic <- var(antithetic_samples)
var_independent <- var(independent_samples)

# 计算方差减少的百分比
percent_reduction <- 100 * (var_independent - var_antithetic) / var_independent

# 输出结果
cat("使用对偶变量的样本方差:", var_antithetic, "\n")
cat("独立样本的方差:", var_independent, "\n")
cat("方差减少百分比:", percent_reduction, "%\n")

## -----------------------------------------------------------------------------
# 设置随机种子以确保可重复性
set.seed(123)

# 定义目标函数 g(x)
g <- function(x) {
  (x^2 / sqrt(2 * pi)) * exp(-x^2 / 2)
}

# 定义重要性函数 f1(x): 正态分布 N(2, 1)
f1_density <- function(x) {
  ifelse(x > 1, dnorm(x, mean = 2, sd = 1), 0)
}

# 定义重要性函数 f2(x): 指数分布 Exp(1)
f2_density <- function(x) {
  ifelse(x > 1, dexp(x - 1, rate = 1), 0)
}

# 采样自 f1(x)
f1_sample <- function(n) {
  x <- rnorm(n, mean = 2, sd = 1)
  x[x > 1] # 只保留 x > 1 的部分
}

# 采样自 f2(x)
f2_sample <- function(n) {
  x <- rexp(n, rate = 1) + 1 # 将指数偏移以从 x = 1 开始
  return(x)
}

# 使用重要性采样估计积分
importance_sampling <- function(n, sample_func, f_density) {
  x <- sample_func(n)
  weights <- g(x) / f_density(x)
  estimate <- mean(weights)
  return(list(estimate = estimate, variance = var(weights)))
}

# 样本数量
n <- 10000

# 估计使用 f1 作为重要性函数
result_f1 <- importance_sampling(n, f1_sample, f1_density)
cat("使用 f1 (N(2, 1)) 的估计值:", result_f1$estimate, "\n")
cat("使用 f1 的方差:", result_f1$variance, "\n")

# 估计使用 f2 作为重要性函数
result_f2 <- importance_sampling(n, f2_sample, f2_density)
cat("使用 f2 (Exp(1)) 的估计值:", result_f2$estimate, "\n")
cat("使用 f2 的方差:", result_f2$variance, "\n")

# 方差比较
variance_reduction <- 100 * (result_f1$variance - result_f2$variance) / result_f1$variance
cat("方差减少百分比:", variance_reduction, "%\n")

## -----------------------------------------------------------------------------
# 加载所需的库
library(ggplot2)

# 定义快速排序的函数
quick_sort <- function(x) {
  if (length(x) <= 1) {
    return(x)
  } else {
    pivot <- x[1]
    less <- x[x < pivot]
    equal <- x[x == pivot]
    greater <- x[x > pivot]
    return(c(quick_sort(less), equal, quick_sort(greater)))
  }
}

# Monte Carlo实验
set.seed(123)  # 设置随机种子，确保结果可重复

# 定义不同的n值
n_values <- c(10^4, 2 * 10^4, 4 * 10^4, 6 * 10^4, 8 * 10^4)

# 模拟次数
n_simulations <- 100

# 创建一个空的数据框来存储结果
results <- data.frame(n = integer(), avg_time = numeric())

# 进行Monte Carlo模拟
for (n in n_values) {
  times <- numeric(n_simulations)
  
  # 对每个n进行100次模拟
  for (i in 1:n_simulations) {
    # 生成1到n的随机排列
    random_numbers <- sample(1:n)
    
    # 计算排序的时间
    start_time <- Sys.time()
    quick_sort(random_numbers)  # 应用快速排序
    end_time <- Sys.time()
    
    # 记录排序时间
    times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  # 计算平均时间
  avg_time <- mean(times)
  
  cat("平均100次模拟的计算时间:", avg_time, "\n", "此时的n:", n, "\n")
  
  # 保存结果
  results <- rbind(results, data.frame(n = n, avg_time = avg_time))
}

# 计算 t_n = n * log(n)
results$t_n <- results$n * log(results$n)

# 回归分析：对 a_n 回归 t_n
fit <- lm(avg_time ~ t_n, data = results)

# 打印回归结果
summary(fit)

# 图形化展示结果，绘制散点图和回归线
ggplot(results, aes(x = t_n, y = avg_time)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(title = "Computation Time vs n log(n)",
       x = "n log(n)",
       y = "Average Computation Time (seconds)") +
  theme_minimal()

