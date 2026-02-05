# nolint: start
#' @title Fix Origin Backtest
#' @description Fixed origin backtesting for AY 2010 with increasing lags.
#'
#' This script performs a fixed origin backtest where:
#' - Training data is updated by adding actual accident year data (not predictions)
#' - Test data is always AY 2010
#' - Initially hold out all but first development lag (m=1)
#' - Increasingly unveil more known lags (m=1 to m=9)
#'
#' Output: Tables of optimal K, lambda, and MAPE, plus data for predictplot.Rmd

library(tidyverse)
library(parallel)
library(doParallel)
library(foreach)
library(glmnet)
library(MASS)
library(caret)

# Source required function files
source("./functions.R")

# =============================================================================
# TUNING FUNCTIONS FOR FIXED ORIGIN BACKTEST
# =============================================================================

#' K tuning function for fixed origin backtest
#' @param m Number of completed lags
#' @param data Training data (fully developed triangles)
#' @param n_fold Number of folds for cross-validation
#' @return List with mean and mean+sd vectors of errors by K
k_tune_fix <- function(m, data, n_fold = 5) {
    # Number of cores available
    no_cores <- detectCores() - 1

    # Filter to accident years with at least m lags complete (AY <= 2010-m)
    # For m=1: use AY <= 2009 (9 years of completed data)
    # For m=9: use AY <= 2001 (1 year of completed data)
    data_filtered <- data %>% filter(accident_year <= 2010 - m)

    # Parallel compute 5-fold cross validation
    error <- mclapply(1:n_fold, function(t) {
        data_fold <- fold_gen(data_filtered, n_fold)
        train_cov <- data_fold$train_cov[[t]]
        test_cov <- data_fold$test_cov[[t]]
        test_matrix <- data_fold$test_matrix[[t]]
        pca <- data_fold$pca[[t]]
        beta <- data_fold$beta[[t]]
        beta_glm <- beta_reg(beta, train_cov, m, alpha = 1, lasso = TRUE)

        error_temp <- rep(0, m)
        for (k in 1:m) {
            pred <- pls_triangle(
                0, m, beta_glm, test_cov, test_matrix, pca, k,
                lasso = TRUE
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

#' Lambda tuning function for fixed origin backtest
#' @param m Number of completed lags
#' @param k Number of principal components
#' @param data Training data
#' @param n_fold Number of folds
#' @return List with optimal lambda and MAPE
l_tune_fix <- function(m, k, data, n_fold = 5) {
    data_filtered <- data %>% filter(accident_year <= 2010 - m)
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
# MAIN PROCESSING FUNCTION FOR EACH LAG
# =============================================================================

#' Process a single lag m in the fixed origin backtest
#' @param m Number of completed lags (1 to 9)
#' @param data_imputed Accumulated training data from prior AYs
#' @param data_cov Full dataset with AY 2010
#' @param Nb Number of bootstrap samples
#' @return List of results including K, lambda, predictions, bootstrap samples
process_m_fix <- function(m, data_imputed, data_cov, Nb = 1000) {
    # (1) Tune K and lambda using historical fully developed triangles
    temp <- k_tune_fix(m, data_imputed)
    K <- which.min(temp$mean)
    temp <- l_tune_fix(m, K, data_imputed)
    lam <- temp$lam
    mape <- temp$mape

    # (2) Generate predictions
    # Training: use all data up to AY 2009 (for training PCA and beta models)
    # Test: always AY 2010
    TT <- data_gen(data_imputed, 2009, 2010) # Train on <= 2009
    TT2 <- data_gen(data_cov, 2009, 2010) # Test on 2010

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

    # Export variables to the cluster
    clusterExport(cl,
        varlist = c(
            "beta", "train_cov", "test_cov",
            "test_matrix", "pca", "K", "lam", "m"
        ),
        envir = environment()
    )

    # Parallel bootstrap loop
    results <- foreach(
        t = 1:Nb,
        .packages = c("dplyr", "glmnet"),
        .export = c(
            "beta_reg", "pls_triangle", "beta", "train_cov", "test_cov",
            "test_matrix", "pca", "K", "lam", "m"
        )
    ) %dopar% {
        bootsize <- 500
        bootid <- sample(nrow(beta), size = bootsize, replace = TRUE)
        beta_boot <- beta[bootid, ]
        cov_boot <- train_cov[bootid, ]
        beta_model_boot <- beta_reg(beta_boot, cov_boot, K, alpha = 1, lasso = TRUE)

        A <- pls_triangle(lam, m, beta_model_boot, test_cov, test_matrix, pca, K)

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

    # (4) Update data_imputed by adding actual AY data (not predictions)
    # Add AY 2002 for m=1, AY 2003 for m=2, etc.
    new_ay <- 2001 + m
    data_imputed_new <- data_cov %>%
        filter(accident_year <= new_ay)

    # Return all results
    list(
        K = K,
        lam = lam,
        mape = mape,
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
        eb = eb, # Residual functions from K+1 to 10th PC
        data_imputed = data_imputed_new,
        n_train = nrow(train_cov)
    )
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

# Initialize parallel computing
num_cores <- parallel::detectCores()
cl <- makeCluster(num_cores - 1)
registerDoParallel(cl)

# Export custom functions to cluster
clusterExport(cl,
    varlist = c("k_tune_fix", "l_tune_fix", "data_gen", "beta_reg", "pls_triangle", "fold_gen"),
    envir = environment()
)

# Initialize tracking variables
results_fix <- list()
Nb_fix <- 1000 # Number of bootstrap samples

# Start with only AY <= 2001 (fully developed triangles)
data_imputed_fix <- data_cov %>% filter(accident_year <= 2001)

start_time <- Sys.time()
cat("Starting fixed origin backtest...\n")

# Process each m value (m = 1 to 9)
# m = 1: only lag 0 known for AY 2010, predict lags 1-9
# m = 9: lags 0-8 known for AY 2010, predict lag 9
for (m in 1:9) {
    cat(sprintf("Processing m = %d (AY 2002 with lag %d)...\n", m, m))
    results_fix[[m]] <- process_m_fix(m, data_imputed_fix, data_cov, Nb_fix)
    data_imputed_fix <- results_fix[[m]]$data_imputed
    cat(sprintf(
        "  K = %d, lambda = %.4f, MAPE = %.4f, n_train = %d\n",
        results_fix[[m]]$K, results_fix[[m]]$lam, results_fix[[m]]$mape, results_fix[[m]]$n_train
    ))
}

end_time <- Sys.time()
runtime <- end_time - start_time
cat("\nTotal runtime:", runtime, "\n")

# Stop the cluster
stopCluster(cl)

# =============================================================================
# EXTRACT RESULTS INTO SEPARATE LISTS
# =============================================================================

k_all_fix <- sapply(results_fix, function(x) x$K)
lam_all_fix <- sapply(results_fix, function(x) x$lam)
mape_all_fix <- sapply(results_fix, function(x) x$mape)
n_train_all_fix <- sapply(results_fix, function(x) x$n_train)
pca_all_fix <- lapply(results_fix, function(x) x$pca)
beta_all_fix <- lapply(results_fix, function(x) x$beta)
beta_model_all_fix <- lapply(results_fix, function(x) x$beta_model)
train_cov_all_fix <- lapply(results_fix, function(x) x$train_cov)
test_cov_all_fix <- lapply(results_fix, function(x) x$test_cov)
test_matrix_all_fix <- lapply(results_fix, function(x) x$test_matrix)
pred_all_fix <- lapply(results_fix, function(x) x$pred)
beta_pred_all_fix <- lapply(results_fix, function(x) x$beta_pred)
betapls_pred_all_fix <- lapply(results_fix, function(x) x$betapls_pred)
boot_pred_all_fix <- lapply(results_fix, function(x) x$boot_pred)
boot_beta_pred_all_fix <- lapply(results_fix, function(x) x$boot_beta_pred)
boot_betapls_pred_all_fix <- lapply(results_fix, function(x) x$boot_betapls_pred)
boot_pred_b_all_fix <- lapply(results_fix, function(x) x$boot_pred_b)
eb_all_fix <- lapply(results_fix, function(x) x$eb) # Residual functions

# =============================================================================
# CREATE OPTIMAL PARAMETERS TABLE
# =============================================================================

optim_param_table_formatted_fix <- data.frame(
    `Accident Year (Lag s)` = paste0(2002:2010, " (", 9:1, ")"),
    `Training Curves` = n_train_all_fix,
    `Optimal K(s)` = k_all_fix,
    `Optimal λ(s)` = round(lam_all_fix, 4),
    `MAPE` = round(mape_all_fix, 4),
    check.names = FALSE
)

cat("\nFormatted Table:\n")
print(optim_param_table_formatted_fix)

# =============================================================================
# SAVE RESULTS
# =============================================================================

save(
    k_all_fix, lam_all_fix, mape_all_fix, n_train_all_fix,
    pca_all_fix, beta_all_fix, beta_model_all_fix, train_cov_all_fix, test_cov_all_fix,
    test_matrix_all_fix, pred_all_fix, betapls_pred_all_fix,
    file = "fix_origin_backtest_output.Rdata"
)

save(
    boot_pred_b_all_fix,
    file = "fix_origin_backtest_bootstrap.Rdata"
)

cat("\nResults saved to:\n")
cat("  - fix_origin_backtest_output.Rdata\n")
cat("  - fix_origin_backtest_bootstrap.Rdata\n")

# =============================================================================
# INTERVAL COVERAGE AND SCORE (using shared parameterized functions)
# =============================================================================

source("./evaluation.R")

cat("\nCalculating coverage and interval scores...\n")
coverage_score_results_fix <- calculate_coverage_score(
    company_list = NULL,
    test_cov_list = test_cov_all_fix,
    test_matrix_list = test_matrix_all_fix,
    pred_list = pred_all_fix,
    boot_pred_b_list = boot_pred_b_all_fix
)

cat("\n=============================================================================\n")
cat("Ultimate Coverage (count and proportion)\n")
cat("=============================================================================\n")
print(coverage_score_results_fix$ultimate_coverage)
cat("\nProportion:\n")
print(round(coverage_score_results_fix$ultimate_coverage_pct, 4))

cat("\n=============================================================================\n")
cat("Functional Coverage (count and proportion)\n")
cat("=============================================================================\n")
print(coverage_score_results_fix$functional_coverage)
cat("\nProportion:\n")
print(round(coverage_score_results_fix$functional_coverage_pct, 4))

cat("\n=============================================================================\n")
cat("Ultimate Interval Score (normalized by EP)\n")
cat("=============================================================================\n")
print(round(coverage_score_results_fix$ultimate_score, 6))

cat("\n=============================================================================\n")
cat("Functional Interval Score (normalized by EP)\n")
cat("=============================================================================\n")
print(round(coverage_score_results_fix$functional_score, 6))

# Save coverage and score results
save(
    coverage_score_results_fix,
    file = "fix_origin_backtest_coverage_score.Rdata"
)
cat("\nCoverage/Score results saved to: fix_origin_backtest_coverage_score.Rdata\n")

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

library(gridExtra)
library(patchwork)
library(fdaoutlier)

# =============================================================================
# PLOTTING FUNCTIONS - MOVED TO OTHER FILES
# =============================================================================
# The following plotting functions have been refactored to separate files:
#   - depth_interval_fix() -> prediction_interval.R
#   - predict_plot_fix() -> forecast_plot.R
#   - predict_plot_grid_fix() -> forecast_plot.R
#   - predict_clr_tracking_fix() -> forecast_plot.R
#
# These functions use the global _fix variables created by this script.
# To use them, ensure prediction_interval.R and forecast_plot.R are sourced.

# =============================================================================
# RESIDUAL FUNCTION (eb) ANALYSIS - i.i.d. TESTING
# =============================================================================

#' Extract and analyze residual functions (eb) for i.i.d. testing
#' eb = beta_high_dim %*% t(pca_high_dim) is the residual from K+1 to 10th PC
#' If the model is correct, eb should be i.i.d. across observations

library(tseries) # For Ljung-Box test

#' Get all residual functions (eb) for a specific m
#' @param m Development lag (1-9)
#' @return Matrix of residual functions (n_obs x 10 lags)
get_eb_fix <- function(m) {
    eb_all_fix[[m]]
}

#' Plot residual functions for visual inspection
#' @param m Development lag
#' @param sample_size Number of curves to plot (NULL = all)
#' @return ggplot object
plot_eb_fix <- function(m, sample_size = 100) {
    eb <- get_eb_fix(m)
    n_obs <- nrow(eb)

    # Sample if too many observations
    if (!is.null(sample_size) && n_obs > sample_size) {
        idx <- sample(n_obs, sample_size)
        eb <- eb[idx, ]
    } else {
        idx <- 1:n_obs
    }

    # Create data frame for plotting
    eb_df <- as.data.frame(eb)
    colnames(eb_df) <- paste0("Lag", 0:9)
    eb_df$obs <- 1:nrow(eb_df)

    eb_long <- eb_df %>%
        pivot_longer(-obs, names_to = "lag", values_to = "value") %>%
        mutate(lag = as.integer(gsub("Lag", "", lag)))

    # Plot
    p <- ggplot(eb_long, aes(x = lag, y = value, group = obs)) +
        geom_line(alpha = 0.3, color = "blue") +
        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            title = element_text(size = 16)
        ) +
        labs(
            title = paste0("Residual Functions (eb) for m = ", m, " (K = ", k_all_fix[m], ")"),
            subtitle = paste0("Showing ", nrow(eb), " of ", n_obs, " observations"),
            x = "Development Lag",
            y = "Residual (eb)"
        ) +
        scale_x_continuous(breaks = 0:9)

    p
}

#' Plot histogram of eb values at each lag
#' @param m Development lag
#' @return ggplot object
plot_eb_hist_fix <- function(m) {
    eb <- get_eb_fix(m)

    eb_df <- as.data.frame(eb)
    colnames(eb_df) <- paste0("Lag ", 0:9)

    eb_long <- eb_df %>%
        pivot_longer(everything(), names_to = "lag", values_to = "value")

    p <- ggplot(eb_long, aes(x = value)) +
        geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.7) +
        facet_wrap(~lag, scales = "free") +
        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 10),
            strip.text = element_text(size = 12)
        ) +
        labs(
            title = paste0("Distribution of Residuals at Each Lag (m = ", m, ")"),
            x = "Residual Value",
            y = "Count"
        )

    p
}

#' Test i.i.d. assumption for residual functions
#' Performs: Shapiro-Wilk (normality), Ljung-Box (autocorrelation), Runs test (randomness)
#' @param m Development lag
#' @return List with test results
test_eb_iid_fix <- function(m) {
    eb <- get_eb_fix(m)
    n_obs <- nrow(eb)
    n_lags <- ncol(eb)

    results <- list()

    # 1. Shapiro-Wilk test for normality (at each lag)
    shapiro_pvals <- numeric(n_lags)
    for (j in 1:n_lags) {
        # Sample if n > 5000 (Shapiro-Wilk limit)
        if (n_obs > 5000) {
            sample_idx <- sample(n_obs, 5000)
            test_data <- eb[sample_idx, j]
        } else {
            test_data <- eb[, j]
        }
        shapiro_pvals[j] <- tryCatch(
            shapiro.test(test_data)$p.value,
            error = function(e) NA
        )
    }
    results$shapiro <- data.frame(
        lag = 0:9,
        p_value = shapiro_pvals,
        reject_05 = shapiro_pvals < 0.05
    )

    # 2. Ljung-Box test for autocorrelation (across observations at each lag)
    ljung_pvals <- numeric(n_lags)
    for (j in 1:n_lags) {
        ljung_pvals[j] <- tryCatch(
            Box.test(eb[, j], lag = min(10, n_obs - 1), type = "Ljung-Box")$p.value,
            error = function(e) NA
        )
    }
    results$ljung_box <- data.frame(
        lag = 0:9,
        p_value = ljung_pvals,
        reject_05 = ljung_pvals < 0.05
    )

    # 3. Runs test for randomness (at each lag)
    runs_pvals <- numeric(n_lags)
    for (j in 1:n_lags) {
        x <- eb[, j]
        med <- median(x)
        runs_pvals[j] <- tryCatch(
            {
                # Convert to binary (above/below median)
                binary <- ifelse(x > med, 1, 0)
                # Calculate runs
                runs <- sum(diff(binary) != 0) + 1
                n1 <- sum(binary == 1)
                n0 <- sum(binary == 0)
                # Expected runs and variance under H0
                mu_r <- (2 * n1 * n0) / (n1 + n0) + 1
                var_r <- (2 * n1 * n0 * (2 * n1 * n0 - n1 - n0)) /
                    ((n1 + n0)^2 * (n1 + n0 - 1))
                # Z-statistic
                z <- (runs - mu_r) / sqrt(var_r)
                2 * pnorm(-abs(z)) # Two-sided p-value
            },
            error = function(e) NA
        )
    }
    results$runs_test <- data.frame(
        lag = 0:9,
        p_value = runs_pvals,
        reject_05 = runs_pvals < 0.05
    )

    # 4. Cross-lag correlation test (test if residuals at different lags are independent)
    cor_matrix <- cor(eb)
    # Extract off-diagonal elements
    off_diag <- cor_matrix[lower.tri(cor_matrix)]
    results$cross_lag_cor <- list(
        cor_matrix = cor_matrix,
        avg_abs_cor = mean(abs(off_diag)),
        max_abs_cor = max(abs(off_diag)),
        min_cor = min(off_diag),
        max_cor = max(off_diag)
    )

    # 5. Summary statistics
    results$summary <- data.frame(
        lag = 0:9,
        mean = colMeans(eb),
        sd = apply(eb, 2, sd),
        skewness = apply(eb, 2, function(x) {
            n <- length(x)
            m3 <- sum((x - mean(x))^3) / n
            s3 <- sd(x)^3
            m3 / s3
        }),
        kurtosis = apply(eb, 2, function(x) {
            n <- length(x)
            m4 <- sum((x - mean(x))^4) / n
            s4 <- sd(x)^4
            m4 / s4 - 3 # Excess kurtosis
        })
    )

    results$m <- m
    results$K <- k_all_fix[m]
    results$n_obs <- n_obs

    class(results) <- "eb_iid_test"
    results
}

#' Print method for eb i.i.d. test results
print.eb_iid_test <- function(x, ...) {
    cat("\n=============================================================================\n")
    cat("Residual Function (eb) i.i.d. Test Results\n")
    cat("m =", x$m, ", K =", x$K, ", n_obs =", x$n_obs, "\n")
    cat("=============================================================================\n")

    cat("\n1. NORMALITY (Shapiro-Wilk Test)\n")
    cat("   H0: Data is normally distributed\n")
    print(x$shapiro, row.names = FALSE)
    cat("   Rejected at 0.05:", sum(x$shapiro$reject_05, na.rm = TRUE), "of 10 lags\n")

    cat("\n2. NO AUTOCORRELATION (Ljung-Box Test)\n")
    cat("   H0: No autocorrelation in the data\n")
    print(x$ljung_box, row.names = FALSE)
    cat("   Rejected at 0.05:", sum(x$ljung_box$reject_05, na.rm = TRUE), "of 10 lags\n")

    cat("\n3. RANDOMNESS (Runs Test)\n")
    cat("   H0: Observations are random\n")
    print(x$runs_test, row.names = FALSE)
    cat("   Rejected at 0.05:", sum(x$runs_test$reject_05, na.rm = TRUE), "of 10 lags\n")

    cat("\n4. CROSS-LAG INDEPENDENCE\n")
    cat("   Average absolute correlation:", round(x$cross_lag_cor$avg_abs_cor, 4), "\n")
    cat("   Max absolute correlation:", round(x$cross_lag_cor$max_abs_cor, 4), "\n")

    cat("\n5. SUMMARY STATISTICS\n")
    print(round(x$summary, 4), row.names = FALSE)

    invisible(x)
}

#' Run all tests for all m values and create summary
test_all_eb_iid_fix <- function() {
    all_tests <- list()
    summary_df <- data.frame(
        m = 1:9,
        K = k_all_fix,
        n_obs = numeric(9),
        shapiro_reject = numeric(9),
        ljung_reject = numeric(9),
        runs_reject = numeric(9),
        avg_cross_cor = numeric(9)
    )

    for (m in 1:9) {
        all_tests[[m]] <- test_eb_iid_fix(m)
        summary_df$n_obs[m] <- all_tests[[m]]$n_obs
        summary_df$shapiro_reject[m] <- sum(all_tests[[m]]$shapiro$reject_05, na.rm = TRUE)
        summary_df$ljung_reject[m] <- sum(all_tests[[m]]$ljung_box$reject_05, na.rm = TRUE)
        summary_df$runs_reject[m] <- sum(all_tests[[m]]$runs_test$reject_05, na.rm = TRUE)
        summary_df$avg_cross_cor[m] <- all_tests[[m]]$cross_lag_cor$avg_abs_cor
    }

    cat("\n=============================================================================\n")
    cat("SUMMARY: Residual Function (eb) i.i.d. Tests Across All m Values\n")
    cat("=============================================================================\n")
    cat("Columns show number of lags (out of 10) where H0 was rejected at α=0.05\n\n")
    print(summary_df, row.names = FALSE)

    list(
        tests = all_tests,
        summary = summary_df
    )
}

# Run the i.i.d. tests
cat("\nRunning i.i.d. tests on residual functions (eb)...\n")
eb_test_results_fix <- test_all_eb_iid_fix()

# Save residual analysis results
save(
    eb_all_fix,
    eb_test_results_fix,
    file = "fix_origin_backtest_eb_analysis.Rdata"
)
cat("\nResidual analysis saved to: fix_origin_backtest_eb_analysis.Rdata\n")

# Example plots (commented out by default)
# plot_eb_fix(5)  # Plot residual curves for m=5
# plot_eb_hist_fix(5)  # Histogram of residuals at each lag
# test_eb_iid_fix(5)  # Detailed test for m=5

# nolint: end
