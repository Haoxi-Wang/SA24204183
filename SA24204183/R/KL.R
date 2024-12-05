#' @title KL function for calculating the KL loss for weighted mixture models.
#' @description This function computes the KL divergence for each model, including AIC, BIC, and weighted models, and performs cross-validation. 
#' It returns adjusted KL losses, weights, and model selection results based on AIC, BIC, and cross-validation.
#' @param i The current model index.
#' @param train The training dataset.
#' @param test The test dataset.
#' @param K The number of candidate models to fit (default is 7).
#' @param R The range for model selection based on AIC/BIC (default is 2).
#' @return A list containing the adjusted KL losses, individual model KL losses, AIC and BIC weights, and CV weights for both 5 and 10 folds \code{n}
#' @import mclust
#' @import distr
#' @import Rsolnp
#' @import dplyr
#' @import mvtnorm
#' @import pgmm
#' @import stats
#' @import Rcpp
#' @useDynLib SA24204183
#' @export
KL_for_sales <- function(i, train, test, K = 7, R = 2) {
  data <- train
  len0 <- length(data)
  ynew <- test
  
  eps <- 1e-06
  epsto0 <- function(x) { x[x < eps] <- 0; x }
  
  fit_GMM_para <- function(data, K, max_iter = 100, tol = 1e-4) {
    # 输入参数：
    # data: 一维数据集 (向量)
    # K: 高斯混合模型中的成分数量
    # max_iter: 最大迭代次数（默认1000）
    # tol: 收敛判定的阈值（默认1e-5）
    if (length(data) < K) stop("K cannot be greater than the number of data points")
    if (any(is.na(data))) stop("Data contains missing values")
    if (var(data) == 0) stop("Data variance is zero")
    
    # 初始化参数
    N <- length(data)
    weights <- rep(1/K, K)  # 初始化权重
    means <- sample(data, K) # 随机初始化均值
    variances <- rep(var(data), K)  # 初始化方差
    sds <- rep(sd(data), K)
    
    # 记录旧的对数似然值，用于判断收敛
    log_likelihood_old <- -Inf
    
    # EM 算法
    for (i in 1:max_iter) {
      # E 步骤：计算每个点属于每个成分的概率
      responsibilities <- matrix(0, nrow = N, ncol = K)
      for (k in 1:K) {
        responsibilities[, k] <- weights[k] * dnorm(data, mean = means[k], sd = sqrt(variances[k]))
      }
      # 归一化得到责任度
      responsibilities <- responsibilities / rowSums(responsibilities)
      
      # M 步骤：更新参数
      N_k <- colSums(responsibilities)  # 每个成分的“软”数据点数量
      for (k in 1:K) {
        # 更新均值
        means[k] <- sum(responsibilities[, k] * data) / N_k[k]
        # 更新方差
        variances[k] <- sum(responsibilities[, k] * (data - means[k])^2) / N_k[k]
        # 更新权重
        weights[k] <- N_k[k] / N
      }
      if (any(is.na(responsibilities)) || any(is.infinite(responsibilities))) stop("Invalid responsibilities")
      if (any(is.na(variances)) || any(variances <= 0)) stop("Invalid variances")
      if (any(is.na(weights)) || any(weights <= 0)) stop("Invalid weights")
      
      # 计算对数似然值
      log_likelihood <- sum(log(rowSums(sapply(1:K, function(k) {
        weights[k] * dnorm(data, mean = means[k], sd = sqrt(variances[k]))
      }))))
      
      # 判断收敛
      if (abs(log_likelihood - log_likelihood_old) < tol) {
        break
      }
      log_likelihood_old <- log_likelihood
    }
    
    # 按均值排序
    order_idx <- order(means)
    weights <- weights[order_idx]
    means <- means[order_idx]
    variances <- variances[order_idx]
    sds <- sqrt(variances)
    
    # 计算 AIC 和 BIC
    num_params <- K * 2 + (K - 1)  # 每个成分均值和方差 + K-1 个独立权重
    AIC <- 2 * num_params - 2 * log_likelihood
    BIC <- num_params * log(N) - 2 * log_likelihood
    
    # 返回结果
    list(
      weights = weights,
      means = means,
      variances = variances,
      sds = sds,
      log_likelihood = log_likelihood,
      AIC = AIC,
      BIC = BIC
    )
  }
  
  
  # 重新定义的安全函数，输出二级清单结果
  fit_GMM_para_safe <- function(data, K, max_iter = 100, tol = 1e-4) {
    tryCatch({
      result <- fit_GMM_para(data, K, max_iter, tol) # 调用原始拟合函数
      if (any(is.na(unlist(result))) || is.null(result)) stop("Invalid output in fit_GMM_para")
      return(list(success = TRUE, result = result)) # 返回二级清单
    }, error = function(e) {
      cat("Error in fit_GMM_para for K =", K, ": ", e$message, "\n")
      return(list(success = FALSE, result = NULL)) # 返回错误信息
    })
  }
  
  myUnivarGMM <- function(means, sds, weights) {
    #means: 均值 sds: 方差 #weights: 权重 
    # 验证输入的长度是否一致
    if (length(means) != length(sds) || length(means) != length(weights)) {
      stop("Means, SDs, and weights must have the same length")
    }
    
    # 归一化混合系数，确保它们的和为1
    normalized_weights <- weights / sum(weights)
    
    # 定义一个返回加权概率密度函数 (PDF) 的函数
    mix_pdf <- function(x) {
      pdf_val <- 0
      
      # 计算每个正态分布的加权 PDF
      for (i in seq_along(means)) {
        pdf_val <- pdf_val + normalized_weights[i] * dnorm(x, mean = means[i], sd = sds[i])
      }
      
      return(pdf_val)
    }
    
    # 返回这个混合分布的 PDF 函数
    return(mix_pdf)
  }
  
  # 1. 尝试拟合各个模型
  thetahat_list <- lapply(1:K, function(k) fit_GMM_para_safe(data, k))
  
  # 检查是否有任何模型拟合失败
  if (any(sapply(thetahat_list, function(x) !x$success))) {
    cat("At least one model fitting failed. Exiting KL_for_sales.\n")
    return(NULL)
  }
  
  # 提取成功的模型结果
  models <- lapply(thetahat_list, function(x) x$result)
  
  # 后续逻辑同之前，使用 `models` 代替原来的 `thetahat_list`
  
  # 定义用于拟合模型的辅助函数
  fit_mixture_models_custom <- function(models_list) {
    lapply(models_list, function(model_params) {
      myUnivarGMM(model_params$means, model_params$sds, model_params$weights)
    })
  }
  
  models_custom <- fit_mixture_models_custom(models)
  
  # 计算 AIC 和 BIC
  AIC_val <- sapply(models, function(info) info$AIC)
  BIC_val <- sapply(models, function(info) info$BIC)
  
  # 根据 AIC 和 BIC 选择候选模型
  P_AIC <- which.min(AIC_val)
  P_BIC <- which.min(BIC_val)
  
  candidates_AIC <- max(1, P_AIC - R):min(K, P_AIC + R)
  candidates_BIC <- max(1, P_BIC - R):min(K, P_BIC + R)
  candidates <- unique(c(candidates_AIC, candidates_BIC))
  len1<-length(candidates)
  
  # 生成候选模型密度函数
  filtered_candidate_models <- models_custom[candidates]
  
  # 计算 KL 损失
  KL_for_each_model <- sapply(1:K, function(i) {
    model <- models_custom[[i]]
    mean(-log(model(ynew)))#KL损失的第二项
  })
  
  #继续使用不同的加权模型方法
  compute_weights <- function(delta_vals) { #辅助计算SAIC和SBIC权重的函数
    #delta_vals: 各个模型对应的\Delta AIC或\Delta BIC值
    exp_vals <- exp(-1 / 2 * delta_vals)
    epsto0(exp_vals / sum(exp_vals))
  }
  
  create_weight_vector <- function(K, candidates, value) {#创建指定位置有给定值、其他位置为0的K维向量
    # 检查 candidates 和 value 的长度是否匹配
    if (length(candidates) != length(value)) {
      stop("The length of 'candidates' and 'value' must be the same.")
    }
    
    # 初始化一个长度为 K 的全零向量
    weight <- numeric(K)
    
    # 将指定位置赋值
    weight[candidates] <- value
    
    # 返回结果
    return(weight)
  }
  
  delta_AIC <- AIC_val - min(AIC_val)
  delta_BIC <- BIC_val - min(BIC_val)
  
  #用于计算KL损失的权重(维数<=K)
  temp_AIC_w <- compute_weights(delta_AIC[candidates])
  temp_BIC_w <- compute_weights(delta_BIC[candidates])
  temp_Equal_w <- c(rep(1/len1,len1))
  
  #用于正式分析权重形态的权重(维数=K)
  AIC_w <- create_weight_vector(K, candidates, temp_AIC_w)
  BIC_w <- create_weight_vector(K, candidates, temp_BIC_w)
  AIC_one=ifelse(seq_along(AIC_w) == which.max(AIC_w), 1, 0)
  BIC_one=ifelse(seq_along(BIC_w) == which.max(BIC_w), 1, 0)
  Equal_w <- create_weight_vector(K,candidates = candidates,c(rep(1/len1,len1)))
  
  # 创建加权混合密度函数
  fit_weighted_model <- function(w, models_set,epsi = eps) {
    #w: 给定权重 
    #models_set: 模型集
    #epsi: 误差
    
    # 确认模型数量
    n_models <- length(models_set)
    
    # 遍历每个模型，检查是否权重接近 1
    for (i in seq_len(n_models)) {
      # 创建一个全零的向量，并在第 i 个位置设置为 1
      test_weights <- rep(0, n_models)
      test_weights[i] <- 1
      
      # 如果权重在第 i 个位置接近 1，直接返回对应模型
      if (all(abs(w - test_weights) < epsi)) {
        return((models_set[[i]]))
      }
    }
    
    # 如果没有发现某个模型权重接近 1，则进行加权平均
    # 定义加权混合密度函数
    mixed_density <- function(x) {
      density_value <- 0  # 初始化密度值为0
      # 逐个遍历 models_set，将各模型的密度乘以相应的 AIC_w 权重并相加
      for (i in seq_len(n_models)) {
        # models_set[[i]] 是一个已经定义好的密度函数，直接传入 x
        density_value <- density_value + w[i] * models_set[[i]](x)
      }
      return(density_value)  # 返回加权混合后的密度值
    }
    
    # 返回最终的加权混合密度函数
    return(mixed_density)
  }
  
  filtered_candidate_models <- models_custom[candidates]
  SAICmodel <- fit_weighted_model(temp_AIC_w, filtered_candidate_models)
  SBICmodel <- fit_weighted_model(temp_BIC_w, filtered_candidate_models)
  Equalmodel <- fit_weighted_model(temp_Equal_w,filtered_candidate_models)
  
  #CV criterion
  calculate_CV_KL <- function(J, data, candidate, ynew) {
    # J: 折数
    # data: 训练集
    # candidate: 筛选后的候选模型个数
    # ynew: 测试集
    
    len0=length(data)     # len0: 数据集长度
    len1=length(candidate)# len1: 候选模型集长度
    
    CV_out <- matrix(nrow = len0, ncol = len1)
    
    for (j in 1:J) {
      # 处理测试集索引边界，确保是整数
      start_index <- floor((j - 1) * len0 / J) + 1
      end_index <- min(ceiling(j * len0 / J), len0)  # 防止超出 len0
      id.test <- start_index:end_index
      id.train <- setdiff(1:len0, id.test)
      
      CV.test <- data[id.test]
      CV.train <- data[id.train]
      
      # 安全拟合 len1 个模型
      CV_thetahat_list <- lapply(1:len1, function(k) fit_GMM_para_safe(CV.train, candidate[k]))
      
      # 检查是否有拟合失败的情况
      if (any(sapply(CV_thetahat_list, function(x) !x$success))) {
        cat("Error: One or more CV model fittings returned NULL.\n")
        return(NULL)
      }
      
      # 提取成功的模型参数
      CV_model_params <- lapply(CV_thetahat_list, function(x) x$result)
      
      # 根据模型参数生成候选模型
      CV_candidate_models <- lapply(CV_model_params, function(params) {
        myUnivarGMM(params$means, params$sds, params$weights)
      })
      
      # 填充 CV_out 矩阵
      for (i in 1:length(CV_candidate_models)) {
        CV_out[start_index:end_index, i] <- CV_candidate_models[[i]](CV.test)
      }
    }
    
    # 定义优化目标函数
    CV.weig <- function(w) {
      -sum(log(CV_out %*% w))
    }
    
    # 初始权重
    w1 <- runif(len1, 0, 1)
    w1 <- w1 / sum(w1)
    
    # 约束条件
    eq <- function(w1) sum(w1)
    ret <- solnp(w1, CV.weig, eqfun = eq, eqB = 1, LB = rep(0, len1), UB = rep(1, len1))
    CV.w <- epsto0(ret$pars)
    CV.w <- CV.w / sum(CV.w)
    
    # 检查是否有单一权重接近 1
    eps <- 1e-6
    for (i in seq_len(len1)) {
      test_weights <- rep(0, len1)
      test_weights[i] <- 1
      
      if (all(abs(CV.w - test_weights) < eps)) {
        CVmodel_mixed <- CV_candidate_models[[i]]
        KL_CV <- mean(-log(CVmodel_mixed(ynew)))
        CV.w <- create_weight_vector(K,candidate,CV.w)
        return(list(KL_CV=KL_CV, CV_w = CV.w))
      }
    }
    
    # 如果没有单一模型占主导，执行加权混合
    CVmodel_mixed <- function(x) {
      density_value <- 0
      for (i in seq_len(len1)) {
        density_value <- density_value + CV.w[i] * CV_candidate_models[[i]](x)
      }
      return(density_value)
    }
    
    KL_CV <- mean(-log(CVmodel_mixed(ynew)))
    CV.w <- create_weight_vector(K,candidate,CV.w)
    
    return(list(KL_CV=KL_CV, CV_w = CV.w))  # 返回 KL_CV 和权重 CV_w
  }
  
  # 计算KL损失
  KL_AIC <- mean(-log((models_custom[[P_AIC]])(ynew)))
  KL_BIC <- mean(-log((models_custom[[P_BIC]])(ynew)))
  KL_SAIC <- mean(-log(SAICmodel(ynew)))
  KL_SBIC <- mean(-log(SBICmodel(ynew)))
  KL_Equal <- mean(-log(Equalmodel(ynew)))
  
  # 调用 calculate_CV_KL，解包结果，分别提取KL损失和权重
  cv_result_5 <- calculate_CV_KL(5, data, candidates, ynew)
  KL_CV5 <- cv_result_5$KL_CV
  CV5_w <- cv_result_5$CV_w
  
  cv_result_10 <- calculate_CV_KL(10, data, candidates, ynew)
  KL_CV10 <- cv_result_10$KL_CV
  CV10_w <- cv_result_10$CV_w
  
  # 返回KL计算结果, 并调整KL损失
  KLloss <- c(KL_AIC, KL_BIC, KL_SAIC, KL_SBIC, KL_Equal,KL_CV5,KL_CV10)
  min_KL_loss <- min(c(KLloss,KL_for_each_model))
  #min_KL_loss=0
  KLloss_adjusted <- KLloss - min_KL_loss  # Adjusted KL loss values
  KL_for_each_model_adjusted <- KL_for_each_model-min_KL_loss
  
  #返回最终运算结果，将 AIC_w, BIC_w 和 CV 权重添加到返回列表中
  list(
    KLloss = KLloss_adjusted,
    KL_for_each_model = KL_for_each_model_adjusted,
    AIC_one=AIC_one,
    BIC_one=BIC_one,
    AIC_w = AIC_w,
    BIC_w = BIC_w,
    Equal_w = Equal_w,
    CV5_w = CV5_w,
    CV10_w = CV10_w
  )
}
