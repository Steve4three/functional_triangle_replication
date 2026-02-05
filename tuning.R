# nolint: start
# Assuming lambda is fixed at 0 for now
k_tune <- function(m, data, n_fold = 5) {
  # Number of cores available
  no_cores <- detectCores() - 1

  # Parallel compute 5-fold cross validation
  error <- mclapply(1:n_fold, function(t) {
    data_fold <- fold_gen(
      data %>% filter(accident_year <= 2010 - m),
      n_fold
    )
    train_cov <- data_fold$train_cov[[t]]
    test_cov <- data_fold$test_cov[[t]]
    test_matrix <- data_fold$test_matrix[[t]]
    pca <- data_fold$pca[[t]]
    beta <- data_fold$beta[[t]]
    beta_glm <- beta_reg(beta, train_cov, m, alpha = 1, lasso = TRUE)

    error_temp <- rep(0, m)
    for (k in 1:m) {
      pred <- pls_triangle(
        0, m, beta_glm, test_cov, test_matrix, pca, k, lasso = TRUE
      )$pred
      pred_diff <- pred - test_matrix[, (m + 1):10]
      pred_sum <- apply(pred_diff, 1, sum)
      test_sum <- apply(test_matrix, 1, sum)
      error_temp[k] <- sum(abs(pred_sum / test_sum)) / nrow(test_matrix)
    }
    error_temp
  }, mc.cores = no_cores)

  mean_vector <- Reduce("+", error) / n_fold
  squared_diffs <- lapply(error, function(x) (x - mean_vector)^2)
  sd_vector <- sqrt(Reduce("+", squared_diffs) / (n_fold - 1))
  list(mean = mean_vector, msd = mean_vector + sd_vector)
}

# Plug in K and tune lambda
l_tune <- function(m, k, data, n_fold = 5) {
  error <- 0
  data_fold <- fold_gen(data, n_fold)
  beta_model <- NULL
  for (t in 1:n_fold) {
    beta_model[[t]] <- beta_reg(
      data_fold$beta[[t]],
      data_fold$train_cov[[t]],
      k,
      alpha = 1,
      lasso = TRUE
    )
  }
  obj <- function(lambda) {
    for (t in 1:n_fold) {
      pred <- pls_triangle(
        lambda, m, beta_model[[t]], data_fold$test_cov[[t]],
        data_fold$test_matrix[[t]], data_fold$pca[[t]], k
      )$pred
      pred_diff <- pred - data_fold$test_matrix[[t]][, (m + 1):10]
      pred_sum <- apply(pred_diff, 1, sum)
      test_sum <- apply(data_fold$test_matrix[[t]], 1, sum)
      error <- error + sum(abs(pred_sum / test_sum)) /
        nrow(data_fold$test_matrix[[t]])
    }
    error / n_fold
  }
  opt <- optimize(obj, c(1e-6, 100))
  list(lam = opt$minimum, mape = opt$objective)
}
# nolint: end