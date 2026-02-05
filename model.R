# nolint: start
#' @title Process m
#' @description Process m
#' @param m Number of completed lags
#' @param data_imputed The imputed data
#' @param data_cov The data
#' @param Nb The number of bootstrap samples
#' @return A list of results
process_m <- function(m, data_imputed, data_cov, Nb = 1000) {
  # (1) Tune K and lambda
  temp <- k_tune(m, data_imputed)
  K <- which.min(temp$mean)
  temp <- l_tune(m, K, data_imputed)
  lam <- temp$lam
  
  # (2) Generate predictions
  TT <- data_gen(data_imputed, 2010-m, 2011-m)
  TT2 <- data_gen(data_cov, 2010-m, 2011-m)
  
  train_cov <- TT$train_cov
  train_matrix <- TT$train_matrix
  test_cov <- TT2$test_cov
  test_matrix <- TT2$test_matrix
  pca <- TT$pca
  beta <- TT$beta
  beta_model <- beta_reg(beta, train_cov, K, alpha = 1, lasso = TRUE)
  
  A <- pls_triangle(lam, m, beta_model, test_cov, test_matrix, pca, K)
  pred <- A$pred
  beta_pred <- A$beta_pred
  betapls_pred <- A$beta_pls
  
  # (3) Bootstrap predictions
  beta_high_dim <- data.matrix(beta[, (K + 1):10])
  pca_high_dim <- data.matrix(pca$rotation[, (K + 1):10])
  eb <- beta_high_dim %*% t(pca_high_dim)
  boot_eb <- array(0, dim = c(nrow(test_matrix), 10 - m, Nb))
  
  for (i in seq_len(nrow(test_matrix))) {
    ebid <- sample(nrow(eb), size = Nb, replace = TRUE)
    boot_eb[i, , ] <- t(eb[ebid, ((m + 1):10)])
  }
  
  # Initialize arrays
  boot_pred <- array(0, dim = c(nrow(test_matrix), 10 - m, Nb))
  boot_beta_pred <- array(0, dim = c(nrow(test_matrix), K, Nb))
  boot_betapls_pred <- array(0, dim = c(nrow(test_matrix), K, Nb))
  
  # Parallel bootstrap
    # Export updated variables to the cluster
  clusterExport(cl, 
    varlist = c("beta", "train_cov", "test_cov", 
    "test_matrix", "pca", "K", "lam", "m"),
    envir = environment()
  )

  # Set up parallel computing for the bootstrap loop 
  results <- foreach(t = 1:Nb, .packages = c("dplyr", "glmnet"), 
    .export = c("beta_reg", "pls_triangle", "beta", "train_cov", "test_cov", 
                "test_matrix", "pca", "K", "lam", "m")) %dopar% {
    bootsize <- 500
    bootid <- sample(nrow(beta), size = bootsize, replace = TRUE)
    beta_boot <- beta[bootid, ] # sample with replacement of beta
    cov_boot <- train_cov[bootid, ] # sample with replacement of covariates
    beta_model_boot <- beta_reg(beta_boot, cov_boot, K, alpha = 1, lasso = TRUE)

    A <- pls_triangle(lam, m, beta_model_boot, test_cov, test_matrix, pca, K)

    list(boot_pred = A$pred,
         boot_beta_pred = A$beta_pred,
         boot_betapls_pred = A$beta_pls)
  }

  # Collect results
  for (t in 1:Nb) {
    boot_pred[, , t] <- results[[t]]$boot_pred
    boot_beta_pred[, , t] <- results[[t]]$boot_beta_pred
    boot_betapls_pred[, , t] <- results[[t]]$boot_betapls_pred
  }
  
  boot_pred_b <- boot_pred + boot_eb
  
  # (4) Update data_imputed
  data_temp <- data_cov %>% filter(accident_year == 2011 - m)
  data_temp[, (8 + m):17] <- pred
  data_imputed <- data_imputed %>% bind_rows(data_temp)
  
  # Return all results
  list(
    K = K,
    lam = lam,
    mape = temp$mape,
    pca = pca,
    beta = beta,
    beta_model = beta_model,
    train_cov = train_cov,
    test_cov = test_cov,
    test_matrix = test_matrix,
    pred = pred,
    beta_pred = beta_pred,
    betapls_pred = betapls_pred,
    boot_pred = boot_pred,
    boot_beta_pred = boot_beta_pred,
    boot_betapls_pred = boot_betapls_pred,
    boot_pred_b = boot_pred_b,
    data_imputed = data_imputed
  )
}


# Initialize parallel computing
num_cores <- parallel::detectCores()
cl <- makeCluster(num_cores - 1)  # Leave one core free
registerDoParallel(cl)

# Alternatively, if your custom functions are in the global environment,
# you can export them
clusterExport(cl, 
  varlist = c("k_tune", "l_tune", "data_gen", "beta_reg", "pls_triangle"),
  envir = environment()
)

# Initialize lists
results <- list()
data_imputed <- data_cov %>% filter(accident_year <= 2001)

start_time <- Sys.time()
# Process each m value
for (m in 9:1) {
  results[[m]] <- process_m(m, data_imputed, data_cov)
  data_imputed <- results[[m]]$data_imputed
}
end_time <- Sys.time()
runtime <- end_time - start_time
cat("\nTotal runtime:", runtime, "\n")

# Stop the cluster when done
stopCluster(cl)

# Extract results into separate lists
k_all <- sapply(results, function(x) x$K)
lam_all <- sapply(results, function(x) x$lam)
mape_all <- sapply(results, function(x) x$mape)
pca_all <- lapply(results, function(x) x$pca)
beta_all <- lapply(results, function(x) x$beta)
beta_model_all <- lapply(results, function(x) x$beta_model)
train_cov_all <- lapply(results, function(x) x$train_cov)
test_cov_all <- lapply(results, function(x) x$test_cov)
test_matrix_all <- lapply(results, function(x) x$test_matrix)
pred_all <- lapply(results, function(x) x$pred)
beta_pred_all <- lapply(results, function(x) x$beta_pred)
betapls_pred_all <- lapply(results, function(x) x$betapls_pred)
boot_pred_all <- lapply(results, function(x) x$boot_pred)
boot_beta_pred_all <- lapply(results, function(x) x$boot_beta_pred)
boot_betapls_pred_all <- lapply(results, function(x) x$boot_betapls_pred)
boot_pred_b_all <- lapply(results, function(x) x$boot_pred_b)

# save(
#   k_all, lam_all, mape_all,
#   pca_all, beta_all, beta_model_all, train_cov_all, test_cov_all,
#   test_matrix_all, pred_all, betapls_pred_all,
#   file = "model_output_out.Rdata"
# )

# save(
#   boot_pred_b_all,
#   file = "model_bootstrap_out.Rdata"
# )

# nolint: end