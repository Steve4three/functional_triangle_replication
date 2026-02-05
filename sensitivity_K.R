# nolint: start
#' @title Sensitivity Analysis for K
#' @description Fixed origin backtesting with fixed K values (2, 4, 6, 8)
#'
#' This script performs sensitivity analysis by:
#' - Fixing K at 2, 4, 6, and 8 (instead of tuning)
#' - Tuning only lambda for each fixed K
#' - Running estimation, bootstrap, coverage and interval score
#' - Comparing with optimal K results from fix_origin_backtest.R
#'
#' Output: Tables comparing MAPE, coverage, and interval scores across K values

library(tidyverse)
library(parallel)
library(doParallel)
library(foreach)
library(glmnet)
library(MASS)
library(caret)
library(fdaoutlier)

# Source required function files
source("./functions.R")

# =============================================================================
# TUNING FUNCTION FOR FIXED K (only tune lambda)
# =============================================================================

#' Lambda tuning function for fixed K
#' @param m Number of completed lags
#' @param k Fixed number of principal components
#' @param data Training data
#' @param n_fold Number of folds
#' @return List with optimal lambda and MAPE
l_tune_sens <- function(m, k, data, n_fold = 5) {
    data_filtered <- data %>% filter(accident_year <= 2010 - m)

    # Check if k is valid (not greater than m)
    if (k > m) {
        return(list(lam = NA, mape = NA))
    }

    data_fold <- fold_gen(data_filtered, n_fold)
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
        error <- 0
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

# =============================================================================
# MAIN PROCESSING FUNCTION FOR FIXED K
# =============================================================================

#' Process a single lag m with fixed K
#' @param m Number of completed lags (1 to 9)
#' @param k Fixed number of principal components
#' @param data_imputed Accumulated training data
#' @param data_cov Full dataset with AY 2010
#' @param Nb Number of bootstrap samples
#' @return List of results
process_m_sens <- function(m, k, data_imputed, data_cov, Nb = 1000) {
    # Check if k is valid
    if (k > m) {
        return(list(
            K = k, lam = NA, mape = NA, valid = FALSE,
            data_imputed = data_cov %>% filter(accident_year <= 2001 + m)
        ))
    }

    # Tune lambda only (K is fixed)
    temp <- l_tune_sens(m, k, data_imputed)
    lam <- temp$lam
    mape <- temp$mape

    # Generate predictions
    TT <- data_gen(data_imputed, 2009, 2010)
    TT2 <- data_gen(data_cov, 2009, 2010)

    train_cov <- TT$train_cov
    train_matrix <- TT$train_matrix
    test_cov <- TT2$test_cov
    test_matrix <- TT2$test_matrix
    pca <- TT$pca
    beta <- TT$beta
    beta_model <- beta_reg(beta, train_cov, k, alpha = 1, lasso = TRUE)

    A <- pls_triangle(lam, m, beta_model, test_cov, test_matrix, pca, k)
    pred <- A$pred
    beta_pred <- A$beta_pred
    betapls_pred <- A$beta_pls

    # Bootstrap predictions
    beta_high_dim <- data.matrix(beta[, (k + 1):10])
    pca_high_dim <- data.matrix(pca$rotation[, (k + 1):10])
    eb <- beta_high_dim %*% t(pca_high_dim)
    boot_eb <- array(0, dim = c(nrow(test_matrix), 10 - m, Nb))

    for (i in seq_len(nrow(test_matrix))) {
        ebid <- sample(nrow(eb), size = Nb, replace = TRUE)
        boot_eb[i, , ] <- t(eb[ebid, ((m + 1):10)])
    }

    # Initialize arrays
    boot_pred <- array(0, dim = c(nrow(test_matrix), 10 - m, Nb))
    boot_beta_pred <- array(0, dim = c(nrow(test_matrix), k, Nb))
    boot_betapls_pred <- array(0, dim = c(nrow(test_matrix), k, Nb))

    # Export variables to the cluster
    clusterExport(cl,
        varlist = c(
            "beta", "train_cov", "test_cov",
            "test_matrix", "pca", "k", "lam", "m"
        ),
        envir = environment()
    )

    # Parallel bootstrap loop
    results <- foreach(
        t = 1:Nb,
        .packages = c("dplyr", "glmnet"),
        .export = c(
            "beta_reg", "pls_triangle", "beta", "train_cov", "test_cov",
            "test_matrix", "pca", "k", "lam", "m"
        )
    ) %dopar% {
        bootsize <- 500
        bootid <- sample(nrow(beta), size = bootsize, replace = TRUE)
        beta_boot <- beta[bootid, ]
        cov_boot <- train_cov[bootid, ]
        beta_model_boot <- beta_reg(beta_boot, cov_boot, k, alpha = 1, lasso = TRUE)

        A <- pls_triangle(lam, m, beta_model_boot, test_cov, test_matrix, pca, k)

        list(
            boot_pred = A$pred,
            boot_beta_pred = A$beta_pred,
            boot_betapls_pred = A$beta_pls
        )
    }

    # Collect results
    for (t in 1:Nb) {
        boot_pred[, , t] <- results[[t]]$boot_pred
        boot_beta_pred[, , t] <- results[[t]]$boot_beta_pred
        boot_betapls_pred[, , t] <- results[[t]]$boot_betapls_pred
    }

    boot_pred_b <- boot_pred + boot_eb

    # Update data_imputed
    new_ay <- 2001 + m
    data_imputed_new <- data_cov %>%
        filter(accident_year <= new_ay)

    list(
        K = k,
        lam = lam,
        mape = mape,
        valid = TRUE,
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
        data_imputed = data_imputed_new,
        n_train = nrow(train_cov)
    )
}

# =============================================================================
# COVERAGE AND SCORE CALCULATION (using shared parameterized functions)
# =============================================================================

source("./evaluation.R")

# =============================================================================
# MAIN EXECUTION - RUN FOR EACH FIXED K VALUE
# =============================================================================

# Fixed K values to test
K_values <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)

# Initialize parallel computing
num_cores <- parallel::detectCores()
cl <- makeCluster(num_cores - 1)
registerDoParallel(cl)

# Export custom functions to cluster
clusterExport(cl,
    varlist = c("l_tune_sens", "data_gen", "beta_reg", "pls_triangle", "fold_gen"),
    envir = environment()
)

# Store results for each K
all_results_sens <- list()
all_coverage_sens <- list()
Nb_sens <- 1000

start_time <- Sys.time()
cat("Starting sensitivity analysis for K...\n")

for (K_fixed in K_values) {
    cat(sprintf("\n=== Processing K = %d ===\n", K_fixed))

    # Initialize for this K
    results_sens <- list()
    data_imputed_sens <- data_cov %>% filter(accident_year <= 2001)

    # Process each m value
    for (m in 1:9) {
        if (K_fixed > m) {
            cat(sprintf("  m = %d: K = %d > m, skipping\n", m, K_fixed))
            results_sens[[m]] <- list(valid = FALSE, K = K_fixed)
            data_imputed_sens <- data_cov %>% filter(accident_year <= 2001 + m)
        } else {
            cat(sprintf("  Processing m = %d with K = %d...\n", m, K_fixed))
            results_sens[[m]] <- process_m_sens(m, K_fixed, data_imputed_sens, data_cov, Nb_sens)
            data_imputed_sens <- results_sens[[m]]$data_imputed
            if (results_sens[[m]]$valid) {
                cat(sprintf(
                    "    lambda = %.4f, MAPE = %.4f\n",
                    results_sens[[m]]$lam, results_sens[[m]]$mape
                ))
            }
        }
    }

    # Extract results for this K
    valid_m <- sapply(results_sens, function(x) x$valid)
    lam_sens <- sapply(results_sens, function(x) if (x$valid) x$lam else NA)
    mape_sens <- sapply(results_sens, function(x) if (x$valid) x$mape else NA)
    n_train_sens <- sapply(results_sens, function(x) if (x$valid) x$n_train else NA)

    test_cov_sens <- lapply(results_sens, function(x) if (x$valid) x$test_cov else NULL)
    test_matrix_sens <- lapply(results_sens, function(x) if (x$valid) x$test_matrix else NULL)
    pred_sens <- lapply(results_sens, function(x) if (x$valid) x$pred else NULL)
    boot_pred_b_sens <- lapply(results_sens, function(x) if (x$valid) x$boot_pred_b else NULL)

    # Calculate coverage and scores using shared function
    coverage_sens <- calculate_coverage_score(
        company_list = NULL,
        test_cov_list = test_cov_sens,
        test_matrix_list = test_matrix_sens,
        pred_list = pred_sens,
        boot_pred_b_list = boot_pred_b_sens,
        valid_m = valid_m
    )

    # Store all results for this K
    all_results_sens[[as.character(K_fixed)]] <- list(
        K = K_fixed,
        valid_m = valid_m,
        lam = lam_sens,
        mape = mape_sens,
        n_train = n_train_sens,
        coverage = coverage_sens
    )

    cat(sprintf("\nResults for K = %d:\n", K_fixed))
    cat("  MAPE: ", paste(round(mape_sens, 4), collapse = ", "), "\n")
    cat("  Coverage (PLS): ", paste(round(coverage_sens$ultimate_coverage_pct[, "PLS"], 4), collapse = ", "), "\n")
}

end_time <- Sys.time()
runtime_sens <- end_time - start_time
cat("\nTotal runtime:", runtime_sens, "\n")

# Stop the cluster
stopCluster(cl)

# =============================================================================
# SAVE SENSITIVITY RESULTS
# =============================================================================

save(
    all_results_sens,
    file = "sensitivity_K_results.Rdata"
)

cat("\nSensitivity results saved to: sensitivity_K_results.Rdata\n")

# nolint: end
