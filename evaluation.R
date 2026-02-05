# =============================================================================
# EVALUATION FUNCTIONS (Parameterized)
# =============================================================================
# Consolidated coverage and interval score calculations.
# Functions accept data lists as parameters for flexibility across scripts.

library(fdaoutlier)

# Source prediction interval functions
source("./prediction_interval.R")

# =============================================================================
# INTERVAL SCORE FUNCTIONS
# =============================================================================

#' Calculate interval score for a single interval
#' @param x Vector of length 2 with lower and upper bounds
#' @param y True value
#' @param alpha Significance level (default 0.05 for 95% interval)
#' @return Interval score
calculate_interval_score <- function(x, y, alpha = 0.05) {
  l <- x[1]
  u <- x[2]
  u - l + 2 / alpha * (
    ifelse(y < l, l - y, 0) %>% as.numeric() +
      ifelse(y > u, y - u, 0) %>% as.numeric()
  )
}

#' Calculate interval score for multiple intervals
#' @param x Matrix with lower bounds in col 1, upper bounds in col 2
#' @param y Vector of true values
#' @param alpha Significance level
#' @return Mean interval score
calculate_interval_score_multiple <- function(x, y, alpha = 0.05) {
  l <- x[, 1]
  u <- x[, 2]
  mean(u - l + 2 / alpha * (
    ifelse(y < l, l - y, 0) + ifelse(y > u, y - u, 0)
  ))
}

# =============================================================================
# COVERAGE AND SCORE CALCULATION
# =============================================================================

#' Calculate coverage and interval score for PLS and EXD methods
#' @param company_list Vector of company SNL keys to evaluate
#' @param test_cov_list List of test covariate data frames (indexed by m)
#' @param test_matrix_list List of test matrices (indexed by m)
#' @param pred_list List of prediction matrices (indexed by m)
#' @param boot_pred_b_list List of bootstrap prediction arrays (indexed by m)
#' @param valid_m Optional logical vector indicating which m values are valid
#' @return List with coverage counts, percentages, and interval scores
calculate_coverage_score <- function(company_list = NULL,
                                     test_cov_list,
                                     test_matrix_list,
                                     pred_list,
                                     boot_pred_b_list,
                                     valid_m = rep(TRUE, 9)) {
  # If no company list provided, find common companies across all valid m
  if (is.null(company_list)) {
    company_list <- NULL
    for (m in 1:9) {
      if (valid_m[m] && !is.null(test_cov_list[[m]])) {
        if (is.null(company_list)) {
          company_list <- test_cov_list[[m]]$snl_key
        } else {
          company_list <- intersect(company_list, test_cov_list[[m]]$snl_key)
        }
      }
    }
  }

  # Handle case where no valid companies
  if (is.null(company_list) || length(company_list) == 0) {
    return(list(
      n_companies = 0,
      ultimate_coverage = matrix(NA, 9, 2, dimnames = list(NULL, c("PLS", "EXD"))),
      ultimate_coverage_pct = matrix(NA, 9, 2, dimnames = list(NULL, c("PLS", "EXD"))),
      functional_coverage = matrix(NA, 9, 2, dimnames = list(NULL, c("PLS", "EXD"))),
      functional_coverage_pct = matrix(NA, 9, 2, dimnames = list(NULL, c("PLS", "EXD"))),
      ultimate_score = matrix(NA, 9, 2, dimnames = list(NULL, c("PLS", "EXD"))),
      functional_score = matrix(NA, 9, 2, dimnames = list(NULL, c("PLS", "EXD")))
    ))
  }

  n_companies <- length(company_list)

  # Initialize result matrices (PLS and EXD only)
  ultimate_coverage <- matrix(0, 9, 2)
  colnames(ultimate_coverage) <- c("PLS", "EXD")

  functional_coverage <- matrix(0, 9, 2)
  colnames(functional_coverage) <- c("PLS", "EXD")

  ultimate_score <- matrix(0, 9, 2)
  colnames(ultimate_score) <- c("PLS", "EXD")

  functional_score <- matrix(0, 9, 2)
  colnames(functional_score) <- c("PLS", "EXD")

  for (m in 1:9) {
    if (!valid_m[m]) {
      ultimate_coverage[m, ] <- NA
      functional_coverage[m, ] <- NA
      ultimate_score[m, ] <- NA
      functional_score[m, ] <- NA
      next
    }

    pls_cov_ult <- 0
    exd_cov_ult <- 0
    pls_cov_func <- 0
    exd_cov_func <- 0
    pls_score_ult <- 0
    exd_score_ult <- 0
    pls_score_func <- 0
    exd_score_func <- 0
    ep_sum <- 0

    for (company in company_list) {
      # Get intervals using parameterized functions
      obj_pls <- pls_interval(
        m, company, test_cov_list, test_matrix_list,
        pred_list, boot_pred_b_list
      )

      if (m < 9) {
        obj_exd <- depth_interval(m, company, test_cov_list, test_matrix_list,
          pred_list, boot_pred_b_list,
          method = "extremal"
        )
      } else {
        obj_exd <- obj_pls
      }

      # Get earned premium for weighting
      ep <- test_cov_list[[m]][
        which(test_cov_list[[m]]$snl_key == company), "ep"
      ]
      ep_sum <- ep_sum + ep

      # True values
      true_ult <- obj_pls$test_c[10]
      true_func <- obj_pls$test_c[(m + 1):10]

      # Ultimate coverage (lag 9 CLR)
      if (m < 9) {
        pls_cov_ult <- pls_cov_ult + ifelse(
          true_ult >= obj_pls$bounds_c95[10 - m, 1] &
            true_ult <= obj_pls$bounds_c95[10 - m, 2], 1, 0
        )
        exd_cov_ult <- exd_cov_ult + ifelse(
          true_ult >= obj_exd$bounds_c95[10 - m, 1] &
            true_ult <= obj_exd$bounds_c95[10 - m, 2], 1, 0
        )

        # Ultimate interval score
        pls_score_ult <- pls_score_ult + calculate_interval_score(
          obj_pls$bounds_c95[10 - m, ], true_ult
        )
        exd_score_ult <- exd_score_ult + calculate_interval_score(
          obj_exd$bounds_c95[10 - m, ], true_ult
        )

        # Functional coverage (all lags covered)
        pls_func_in <- sum(true_func < obj_pls$bounds_c95[, 1]) +
          sum(true_func > obj_pls$bounds_c95[, 2])
        exd_func_in <- sum(true_func < obj_exd$bounds_c95[, 1]) +
          sum(true_func > obj_exd$bounds_c95[, 2])

        pls_cov_func <- pls_cov_func + ifelse(pls_func_in == 0, 1, 0)
        exd_cov_func <- exd_cov_func + ifelse(exd_func_in == 0, 1, 0)

        # Functional interval score
        pls_score_func <- pls_score_func + calculate_interval_score_multiple(
          obj_pls$bounds_c95, true_func
        )
        exd_score_func <- exd_score_func + calculate_interval_score_multiple(
          obj_exd$bounds_c95, true_func
        )
      } else {
        # m = 9: single value
        pls_cov_ult <- pls_cov_ult + ifelse(
          true_ult >= obj_pls$bounds_c95[1] &
            true_ult <= obj_pls$bounds_c95[2], 1, 0
        )
        exd_cov_ult <- exd_cov_ult + ifelse(
          true_ult >= obj_exd$bounds_c95[1] &
            true_ult <= obj_exd$bounds_c95[2], 1, 0
        )

        pls_score_ult <- pls_score_ult + calculate_interval_score(
          obj_pls$bounds_c95, true_ult
        )
        exd_score_ult <- exd_score_ult + calculate_interval_score(
          obj_exd$bounds_c95, true_ult
        )

        pls_cov_func <- pls_cov_func + ifelse(
          true_func >= obj_pls$bounds_c95[1] &
            true_func <= obj_pls$bounds_c95[2], 1, 0
        )
        exd_cov_func <- exd_cov_func + ifelse(
          true_func >= obj_exd$bounds_c95[1] &
            true_func <= obj_exd$bounds_c95[2], 1, 0
        )

        pls_score_func <- pls_score_func + calculate_interval_score(
          obj_pls$bounds_c95, true_func
        )
        exd_score_func <- exd_score_func + calculate_interval_score(
          obj_exd$bounds_c95, true_func
        )
      }
    }

    # Store results
    ultimate_coverage[m, ] <- c(pls_cov_ult, exd_cov_ult)
    functional_coverage[m, ] <- c(pls_cov_func, exd_cov_func)
    ultimate_score[m, ] <- c(pls_score_ult, exd_score_ult) / ep_sum
    functional_score[m, ] <- c(pls_score_func, exd_score_func) / ep_sum
  }

  # Convert coverage counts to proportions
  ultimate_coverage_pct <- ultimate_coverage / n_companies
  functional_coverage_pct <- functional_coverage / n_companies

  list(
    n_companies = n_companies,
    ultimate_coverage = ultimate_coverage,
    ultimate_coverage_pct = ultimate_coverage_pct,
    functional_coverage = functional_coverage,
    functional_coverage_pct = functional_coverage_pct,
    ultimate_score = ultimate_score,
    functional_score = functional_score
  )
}
# =============================================================================
# ULTIMATE AND FUNCTIONAL METRICS CALCULATION (PLS, EXD, BD, CL)
# =============================================================================

#' Calculate ultimate and functional coverage metrics for all methods
#' This computes MAPE from actual predictions, not from tuning
#' @param company_list Vector of company SNL keys to evaluate
#' @param test_cov_list List of test covariate data frames
#' @param test_matrix_list List of test matrices
#' @param pred_list List of prediction matrices
#' @param boot_pred_b_list List of bootstrap prediction arrays
#' @return List with detailed tables and summary statistics
calculate_ultimate_metrics <- function(company_list,
                                       test_cov_list,
                                       test_matrix_list,
                                       pred_list,
                                       boot_pred_b_list) {
  # Initialize summary matrices
  # Columns: MAPE (2), Ultimate coverage (4), Ultimate score (4), Functional coverage (4), Functional score (4)
  summary_table <- matrix(0, 9, 18)
  colnames(summary_table) <- c(
    "PLS_MAPE", "CL_MAPE",
    "PLS_cov", "EXD_cov", "BD_cov", "CL_cov",
    "PLS_score", "EXD_score", "BD_score", "CL_score",
    "PLS_func_cov", "EXD_func_cov", "BD_func_cov", "CL_func_cov",
    "PLS_func_score", "EXD_func_score", "BD_func_score", "CL_func_score"
  )

  n_companies <- length(company_list)
  ultimate_tables <- list()

  for (m in 1:9) {
    # Initialize table for this lag
    table_data <- data.frame(
      company = character(),
      true_ult = numeric(),
      pls_pred = numeric(),
      cl_pred = numeric(),
      pls_lb = numeric(),
      pls_ub = numeric(),
      exd_lb = numeric(),
      exd_ub = numeric(),
      bd_lb = numeric(),
      bd_ub = numeric(),
      cl_lb = numeric(),
      cl_ub = numeric(),
      pls_err = numeric(),
      cl_err = numeric(),
      pls_cov = numeric(),
      exd_cov = numeric(),
      bd_cov = numeric(),
      cl_cov = numeric(),
      pls_func_cov = numeric(),
      exd_func_cov = numeric(),
      bd_func_cov = numeric(),
      cl_func_cov = numeric(),
      stringsAsFactors = FALSE
    )

    # Ultimate scores (sum then divide by ep_sum)
    pls_score_sum <- 0
    exd_score_sum <- 0
    bd_score_sum <- 0
    cl_score_sum <- 0
    # Functional scores
    pls_func_score_sum <- 0
    exd_func_score_sum <- 0
    bd_func_score_sum <- 0
    cl_func_score_sum <- 0
    ep_sum <- 0
    i <- 1

    for (company in company_list) {
      # Check if company is in test set for this m
      if (!(company %in% test_cov_list[[m]]$snl_key)) next

      # Get PLS interval (returns values in DOLLARS, scaled by ep)
      obj_pls <- tryCatch(
        pls_interval(
          m, company, test_cov_list, test_matrix_list,
          pred_list, boot_pred_b_list
        ),
        error = function(e) NULL
      )
      if (is.null(obj_pls)) next

      # Get depth intervals (EXD and BD) - also in DOLLARS
      if (m < 9) {
        obj_exd <- tryCatch(
          depth_interval(m, company, test_cov_list, test_matrix_list,
            pred_list, boot_pred_b_list,
            method = "extremal"
          ),
          error = function(e) obj_pls
        )
        obj_bd <- tryCatch(
          depth_interval(m, company, test_cov_list, test_matrix_list,
            pred_list, boot_pred_b_list,
            method = "bd"
          ),
          error = function(e) obj_pls
        )
      } else {
        obj_exd <- obj_pls
        obj_bd <- obj_pls
      }

      # Get CL predictions (returns in DOLLARS)
      obj_cl <- tryCatch(cl(company), error = function(e) NULL)
      if (is.null(obj_cl)) next

      # Get earned premium for weighting
      ep <- test_cov_list[[m]][which(test_cov_list[[m]]$snl_key == company), "ep"]
      ep_sum <- ep_sum + ep

      # True values (in DOLLARS from pls_interval)
      true_ult <- obj_pls$test_c[10]
      true_func <- obj_pls$test_c[(m + 1):10]

      # === ULTIMATE METRICS ===
      if (m < 9) {
        pls_pred <- obj_pls$pred_c[10 - m]
        pls_bounds <- obj_pls$bounds_c95[10 - m, ]
        exd_bounds <- obj_exd$bounds_c95[10 - m, ]
        bd_bounds <- obj_bd$bounds_c95[10 - m, ]
      } else {
        pls_pred <- obj_pls$pred_c[1]
        pls_bounds <- obj_pls$bounds_c95
        exd_bounds <- obj_exd$bounds_c95
        bd_bounds <- obj_bd$bounds_c95
      }

      # CL predicted ultimate (already in DOLLARS - no division by ep!)
      cl_pred <- obj_cl$pred[[m]][10 - m]
      cl_se <- sqrt(obj_cl$msep[[m]][10 - m])
      cl_lb <- cl_pred - qnorm(0.975) * cl_se
      cl_ub <- cl_pred + qnorm(0.975) * cl_se

      # Errors
      pls_err <- pls_pred - true_ult
      cl_err <- cl_pred - true_ult

      # Ultimate coverage indicators
      pls_cov <- ifelse(true_ult >= pls_bounds[1] & true_ult <= pls_bounds[2], 1, 0)
      exd_cov <- ifelse(true_ult >= exd_bounds[1] & true_ult <= exd_bounds[2], 1, 0)
      bd_cov <- ifelse(true_ult >= bd_bounds[1] & true_ult <= bd_bounds[2], 1, 0)
      cl_cov <- ifelse(true_ult >= cl_lb & true_ult <= cl_ub, 1, 0)

      # Ultimate interval scores (sum, will divide by ep_sum later)
      pls_score_sum <- pls_score_sum + calculate_interval_score(pls_bounds, true_ult)
      exd_score_sum <- exd_score_sum + calculate_interval_score(exd_bounds, true_ult)
      bd_score_sum <- bd_score_sum + calculate_interval_score(bd_bounds, true_ult)
      cl_score_sum <- cl_score_sum + calculate_interval_score(c(cl_lb, cl_ub), true_ult)

      # === FUNCTIONAL METRICS ===
      if (m < 9) {
        # PLS/EXD/BD functional bounds (matrix: rows = lags, cols = lower/upper)
        pls_func_bounds <- obj_pls$bounds_c95
        exd_func_bounds <- obj_exd$bounds_c95
        bd_func_bounds <- obj_bd$bounds_c95

        # CL functional bounds (in DOLLARS - no division by ep)
        cl_pred_vec <- obj_cl$pred[[m]] # vector of length 10-m
        cl_se_vec <- sqrt(obj_cl$msep[[m]])
        cl_lb_vec <- cl_pred_vec - qnorm(0.975) * cl_se_vec
        cl_ub_vec <- cl_pred_vec + qnorm(0.975) * cl_se_vec

        # Functional coverage: all points within bounds
        pls_func_in <- sum(true_func < pls_func_bounds[, 1]) +
          sum(true_func > pls_func_bounds[, 2])
        exd_func_in <- sum(true_func < exd_func_bounds[, 1]) +
          sum(true_func > exd_func_bounds[, 2])
        bd_func_in <- sum(true_func < bd_func_bounds[, 1]) +
          sum(true_func > bd_func_bounds[, 2])
        cl_func_in <- sum(true_func < cl_lb_vec) + sum(true_func > cl_ub_vec)

        pls_func_cov <- ifelse(pls_func_in == 0, 1, 0)
        exd_func_cov <- ifelse(exd_func_in == 0, 1, 0)
        bd_func_cov <- ifelse(bd_func_in == 0, 1, 0)
        cl_func_cov <- ifelse(cl_func_in == 0, 1, 0)

        # Functional interval scores (using mean from calculate_interval_score_multiple)
        pls_func_score_sum <- pls_func_score_sum +
          calculate_interval_score_multiple(pls_func_bounds, true_func)
        exd_func_score_sum <- exd_func_score_sum +
          calculate_interval_score_multiple(exd_func_bounds, true_func)
        bd_func_score_sum <- bd_func_score_sum +
          calculate_interval_score_multiple(bd_func_bounds, true_func)
        cl_func_bounds_mat <- cbind(cl_lb_vec, cl_ub_vec)
        cl_func_score_sum <- cl_func_score_sum +
          calculate_interval_score_multiple(cl_func_bounds_mat, true_func)
      } else {
        # m = 9: functional = ultimate (single point)
        pls_func_cov <- pls_cov
        exd_func_cov <- exd_cov
        bd_func_cov <- bd_cov
        cl_func_cov <- cl_cov

        pls_func_score_sum <- pls_func_score_sum +
          calculate_interval_score(pls_bounds, true_func)
        exd_func_score_sum <- exd_func_score_sum +
          calculate_interval_score(exd_bounds, true_func)
        bd_func_score_sum <- bd_func_score_sum +
          calculate_interval_score(bd_bounds, true_func)
        cl_func_score_sum <- cl_func_score_sum +
          calculate_interval_score(c(cl_lb, cl_ub), true_func)
      }

      # Store row
      table_data[i, ] <- list(
        company, true_ult, pls_pred, cl_pred,
        pls_bounds[1], pls_bounds[2],
        exd_bounds[1], exd_bounds[2],
        bd_bounds[1], bd_bounds[2],
        cl_lb, cl_ub,
        pls_err, cl_err,
        pls_cov, exd_cov, bd_cov, cl_cov,
        pls_func_cov, exd_func_cov, bd_func_cov, cl_func_cov
      )
      i <- i + 1
    }

    ultimate_tables[[m]] <- table_data

    # Calculate MAPE (percentage)
    if (nrow(table_data) > 0) {
      summary_table[m, "PLS_MAPE"] <- mean(abs(table_data$pls_err / table_data$true_ult), na.rm = TRUE)
      summary_table[m, "CL_MAPE"] <- mean(abs(table_data$cl_err / table_data$true_ult), na.rm = TRUE)

      # Ultimate coverage percentages
      summary_table[m, "PLS_cov"] <- sum(table_data$pls_cov, na.rm = TRUE) / n_companies
      summary_table[m, "EXD_cov"] <- sum(table_data$exd_cov, na.rm = TRUE) / n_companies
      summary_table[m, "BD_cov"] <- sum(table_data$bd_cov, na.rm = TRUE) / n_companies
      summary_table[m, "CL_cov"] <- sum(table_data$cl_cov, na.rm = TRUE) / n_companies

      # Ultimate interval scores (sum / ep_sum as per old code)
      summary_table[m, "PLS_score"] <- pls_score_sum / ep_sum
      summary_table[m, "EXD_score"] <- exd_score_sum / ep_sum
      summary_table[m, "BD_score"] <- bd_score_sum / ep_sum
      summary_table[m, "CL_score"] <- cl_score_sum / ep_sum

      # Functional coverage percentages
      summary_table[m, "PLS_func_cov"] <- sum(table_data$pls_func_cov, na.rm = TRUE) / n_companies
      summary_table[m, "EXD_func_cov"] <- sum(table_data$exd_func_cov, na.rm = TRUE) / n_companies
      summary_table[m, "BD_func_cov"] <- sum(table_data$bd_func_cov, na.rm = TRUE) / n_companies
      summary_table[m, "CL_func_cov"] <- sum(table_data$cl_func_cov, na.rm = TRUE) / n_companies

      # Functional interval scores (sum / ep_sum)
      summary_table[m, "PLS_func_score"] <- pls_func_score_sum / ep_sum
      summary_table[m, "EXD_func_score"] <- exd_func_score_sum / ep_sum
      summary_table[m, "BD_func_score"] <- bd_func_score_sum / ep_sum
      summary_table[m, "CL_func_score"] <- cl_func_score_sum / ep_sum
    }
  }

  list(
    tables = ultimate_tables,
    summary = summary_table,
    n_companies = n_companies
  )
}

# =============================================================================
# PIT PLOT FUNCTIONS
# =============================================================================

#' PIT histogram plot for a single m and lag (parameterized version)
#' @param m Number of completed lags
#' @param j Lag index (1 to 10-m)
#' @param company_list Vector of company SNL keys
#' @param test_cov_list List of test covariate data frames
#' @param test_matrix_list List of test matrices
#' @param boot_pred_b_list List of bootstrap prediction arrays
#' @return List with plot and KS statistic
pit_plot <- function(m, j, company_list, test_cov_list, test_matrix_list,
                     boot_pred_b_list) {
  ind <- which(test_cov_list[[m]]$snl_key %in% company_list)
  test_matrix <- test_matrix_list[[m]][ind, ]
  test_matrix_c <- t(apply(test_matrix[, 1:10], 1, cumsum))
  quants1 <- array(0, dim = c(nrow(test_matrix), 10 - m))
  boot_pred_b <- boot_pred_b_list[[m]][ind, , ]

  for (n in seq_len(nrow(test_matrix))) {
    if (m < 9) {
      zz <- test_matrix_c[n, m] + t(apply(boot_pred_b[n, , ], 2, cumsum))
    } else {
      zz <- test_matrix_c[n, m] + t(boot_pred_b[n, , ])
    }
    for (s in (m + 1):10) {
      quants1[n, s - m] <- sum(zz[, s - m] < test_matrix_c[n, s]) / 1000
    }
  }

  plot_data <- data.frame(x = quants1[, j])
  x_label <- paste0(
    "Predictive quantiles of CLR at Lag ", m + j - 1,
    " for ", m, " complete lags"
  )

  p1 <- ggplot(data = plot_data, aes(x = x)) +
    geom_histogram(bins = 10) +
    theme_bw() +
    xlab(x_label) +
    ylab("Frequency")

  sorted_quants <- sort(quants1[, j])
  n_obs <- nrow(test_matrix)
  theoretical_quantiles <- seq(1 / n_obs, 1, length.out = n_obs)
  ks <- max(abs(sorted_quants - theoretical_quantiles))

  list(plt = p1, ks = ks)
}

#' PIT histogram plot pooling across m=3:6 for lags 7-9 (legacy version)
#' This version is used by main.R for the ECDF calibration plot (Figure 7)
#' Note: This version requires global variables: test_cov_all, test_matrix_all,
#' boot_pred_b_all, and colist (from main.R execution context)
#' @param j Lag index (7, 8, or 9)
#' @return List with quants, plot, KS statistic, and number of observations
pit_plot2 <- function(j) {
  quants <- numeric()
  for (m in 3:6) {
    ind <- which(test_cov_all[[m]]$snl_key %in% colist)
    test_matrix <- test_matrix_all[[m]][ind, ]
    test_matrix_c <- t(apply(test_matrix[, 1:10], 1, cumsum))
    quants1 <- array(0, dim = c(nrow(test_matrix), 10 - m))
    boot_pred_b <- boot_pred_b_all[[m]][ind, , ]
    for (n in seq_len(nrow(test_matrix))) {
      if (m < 9) {
        zz <- test_matrix_c[n, m] + t(apply(boot_pred_b[n, , ], 2, cumsum))
      } else {
        zz <- test_matrix_c[n, m] + t(boot_pred_b[n, , ])
      }
      for (s in (m + 1):10) {
        quants1[n, s - m] <- sum(zz[, s - m] < test_matrix_c[n, s]) / 1000
      }
    }
    quants <- quants %>% rbind(quants1[, (7 - m):(9 - m)])
  }

  plot_data <- data.frame(x = quants[, j - 6])
  x_label <- paste0("Pred. quantiles of CLR at Lag ", j)

  plt <- ggplot(data = plot_data, aes(x = x)) +
    geom_histogram(bins = 10) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text = element_text(size = 24),
      axis.title = element_text(size = 24),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    labs(x = x_label, y = "Frequency")

  sorted_quants <- sort(quants[, j - 6])
  n_obs <- nrow(quants)
  theoretical_quantiles <- seq(1 / n_obs, 1, length.out = n_obs)
  ks <- max(abs(sorted_quants - theoretical_quantiles))

  list(quants = quants[, j - 6], plt = plt, ks = ks, n = nrow(quants))
}

# =============================================================================
# ECDF PLOT HELPER
# =============================================================================

#' Create ECDF plot for calibration assessment
#' @param quants Vector of predictive quantiles
#' @param lag Lag value for x-axis label
#' @param show_y_title Whether to show y-axis title (TRUE for first plot in row)
#' @return ggplot object
create_ecdf_plot <- function(quants, lag, show_y_title = FALSE) {
  ggplot(data.frame(q = quants), aes(x = q)) +
    stat_ecdf(geom = "step", size = 0.8, pad = FALSE) +
    geom_line(
      data = data.frame(x = c(0, 1), y = c(0, 1)),
      aes(x = x, y = y), color = "red", linetype = "dashed", size = 0.8
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text = element_text(size = 24),
      axis.title.x = element_text(size = 24),
      axis.title.y = if (show_y_title) element_text(size = 24) else element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    labs(
      x = paste0("Pred. quantiles of Lag ", lag, " CLR"),
      y = if (show_y_title) "Cumulative Probability" else NULL
    )
}
