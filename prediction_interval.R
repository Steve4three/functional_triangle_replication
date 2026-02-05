# =============================================================================
# PREDICTION INTERVAL FUNCTIONS (Parameterized)
# =============================================================================
# These functions accept data lists as parameters rather than relying on
# global variables, making them reusable across different backtesting scripts.

library(fdaoutlier)

#' Pointwise prediction interval
#' @param m Number of completed lags
#' @param company Company SNL key
#' @param test_cov_list List of test covariate data frames (indexed by m)
#' @param test_matrix_list List of test matrices (indexed by m)
#' @param pred_list List of prediction matrices (indexed by m)
#' @param boot_pred_b_list List of bootstrap prediction arrays (indexed by m)
#' @return List with test, prediction, and bounds data (scaled by ep)
pls_interval <- function(m, company, test_cov_list, test_matrix_list,
                         pred_list, boot_pred_b_list) {
  n <- which(test_cov_list[[m]]$snl_key == company)
  ep <- test_cov_list[[m]][n, "ep"]
  test <- test_matrix_list[[m]][n, ]
  test_c <- cumsum(test)

  A <- boot_pred_b_list[[m]][n, , ]
  pred <- pred_list[[m]][n, ]
  pred_c <- test_c[m] + cumsum(pred)

  if (m < 9) {
    boot_i <- t(A)
    boot_c <- apply(boot_i, 1, cumsum) %>% t()

    probs <- c(0.25, 0.75, 0.05, 0.95, 0.025, 0.975)
    bounds_i <- apply(boot_i, 2, quantile, probs = probs)
    bounds_c <- test_c[m] + apply(bounds_i, 1, cumsum)

    bounds_i50 <- bounds_i[1:2, ] %>% t()
    bounds_i90 <- bounds_i[3:4, ] %>% t()
    bounds_i95 <- bounds_i[5:6, ] %>% t()
    bounds_c50 <- bounds_c[, 1:2]
    bounds_c90 <- bounds_c[, 3:4]
    bounds_c95 <- bounds_c[, 5:6]
  } else {
    boot_i <- A
    boot_c <- boot_i

    probs <- c(0.25, 0.75, 0.05, 0.95, 0.025, 0.975)
    bounds_i <- quantile(A, probs = probs)
    bounds_c <- test_c[m] + bounds_i

    bounds_i50 <- bounds_i[1:2]
    bounds_i90 <- bounds_i[3:4]
    bounds_i95 <- bounds_i[5:6]
    bounds_c50 <- bounds_c[1:2]
    bounds_c90 <- bounds_c[3:4]
    bounds_c95 <- bounds_c[5:6]
  }

  list(
    test_i = test * ep, test_c = test_c * ep,
    pred_i = pred * ep, pred_c = pred_c * ep,
    boot_i = boot_i * ep, boot_c = boot_c * ep,
    bounds_i50 = bounds_i50 * ep, bounds_i90 = bounds_i90 * ep,
    bounds_i95 = bounds_i95 * ep, bounds_c50 = bounds_c50 * ep,
    bounds_c90 = bounds_c90 * ep, bounds_c95 = bounds_c95 * ep
  )
}

#' Depth-based prediction interval
#' @param m Number of completed lags
#' @param company Company SNL key
#' @param test_cov_list List of test covariate data frames (indexed by m)
#' @param test_matrix_list List of test matrices (indexed by m)
#' @param pred_list List of prediction matrices (indexed by m)
#' @param boot_pred_b_list List of bootstrap prediction arrays (indexed by m)
#' @param method Depth method: "extremal", "mbd", "tvd", "bd"
#' @param ef Empirical factor for functional boxplot
#' @return List with test, prediction, and bounds data (scaled by ep)
depth_interval <- function(m, company, test_cov_list, test_matrix_list,
                           pred_list, boot_pred_b_list,
                           method = "extremal", ef = 2.5) {
  n <- which(test_cov_list[[m]]$snl_key == company)
  ep <- test_cov_list[[m]][n, "ep"]
  test <- test_matrix_list[[m]][n, ]
  test_c <- cumsum(test)

  A <- boot_pred_b_list[[m]][n, , ]
  pred <- pred_list[[m]][n, ]
  pred_c <- test_c[m] + cumsum(pred)

  if (m < 9) {
    boot_i <- t(A)
    boot_c <- apply(boot_i, 1, cumsum) %>% t()
    depth_val <- functional_boxplot(
      boot_i,
      depth_method = method, emp_factor = ef
    )$depth_values
    curve50 <- boot_i[depth_val > sort(depth_val)[dim(boot_i)[1] * 0.5], ]
    curve90 <- boot_i[depth_val > sort(depth_val)[dim(boot_i)[1] * 0.1], ]
    curve95 <- boot_i[depth_val > sort(depth_val)[dim(boot_i)[1] * 0.05], ]
    bounds_i50 <- cbind(apply(curve50, 2, min), apply(curve50, 2, max))
    bounds_i90 <- cbind(apply(curve90, 2, min), apply(curve90, 2, max))
    bounds_i95 <- cbind(apply(curve95, 2, min), apply(curve95, 2, max))
    depth_val <- functional_boxplot(
      boot_c,
      depth_method = method, emp_factor = ef
    )$depth_values
    curve50 <- boot_c[depth_val > sort(depth_val)[dim(boot_c)[1] * 0.5], ]
    curve90 <- boot_c[depth_val > sort(depth_val)[dim(boot_c)[1] * 0.1], ]
    curve95 <- boot_c[depth_val > sort(depth_val)[dim(boot_c)[1] * 0.05], ]
    bounds_c50 <- test_c[m] + cbind(
      apply(curve50, 2, min), apply(curve50, 2, max)
    )
    bounds_c90 <- test_c[m] + cbind(
      apply(curve90, 2, min), apply(curve90, 2, max)
    )
    bounds_c95 <- test_c[m] + cbind(
      apply(curve95, 2, min), apply(curve95, 2, max)
    )
  } else {
    # m = 9 use pointwise quantile for the last lag
    boot_i <- A
    boot_c <- boot_i
    bounds_i50 <- quantile(A, probs = c(0.25, 0.75))
    bounds_i90 <- quantile(A, probs = c(0.05, 0.95))
    bounds_i95 <- quantile(A, probs = c(0.025, 0.975))
    bounds_c50 <- test_c[m] + bounds_i50
    bounds_c90 <- test_c[m] + bounds_i90
    bounds_c95 <- test_c[m] + bounds_i95
  }

  list(
    test_i = test * ep, test_c = test_c * ep,
    pred_i = pred * ep, pred_c = pred_c * ep,
    boot_i = boot_i * ep, boot_c = boot_c * ep,
    bounds_i50 = bounds_i50 * ep, bounds_i90 = bounds_i90 * ep,
    bounds_i95 = bounds_i95 * ep, bounds_c50 = bounds_c50 * ep,
    bounds_c90 = bounds_c90 * ep, bounds_c95 = bounds_c95 * ep
  )
}

# =============================================================================
# FIXED ORIGIN BACKTEST INTERVAL FUNCTION
# =============================================================================
# Uses global variables: test_cov_all_fix, test_matrix_all_fix,
# pred_all_fix, boot_pred_b_all_fix

#' Depth-based prediction interval for fixed origin backtest
#' Uses the _fix suffix global variables
#' @param m Number of completed lags
#' @param company Company SNL key
#' @param method Depth method: "extremal", "mbd", "tvd", "bd"
#' @param ef Empirical factor for functional boxplot
#' @return List with test, prediction, and bounds data
depth_interval_fix <- function(m, company, method = "extremal", ef = 2.5) {
  n <- which(test_cov_all_fix[[m]]$snl_key == company)
  ep <- test_cov_all_fix[[m]][n, "ep"]
  test <- test_matrix_all_fix[[m]][n, ]
  test_c <- cumsum(test)

  A <- boot_pred_b_all_fix[[m]][n, , ]
  pred <- pred_all_fix[[m]][n, ]
  pred_c <- test_c[m] + cumsum(pred)

  if (m < 9) {
    boot_i <- t(A)
    boot_c <- apply(boot_i, 1, cumsum) %>% t()
    depth_val <- functional_boxplot(
      boot_i,
      depth_method = method, emp_factor = ef
    )$depth_values
    curve50 <- boot_i[depth_val > sort(depth_val)[dim(boot_i)[1] * 0.5], ]
    curve90 <- boot_i[depth_val > sort(depth_val)[dim(boot_i)[1] * 0.1], ]
    curve95 <- boot_i[depth_val > sort(depth_val)[dim(boot_i)[1] * 0.05], ]
    bounds_i50 <- cbind(apply(curve50, 2, min), apply(curve50, 2, max))
    bounds_i90 <- cbind(apply(curve90, 2, min), apply(curve90, 2, max))
    bounds_i95 <- cbind(apply(curve95, 2, min), apply(curve95, 2, max))
    depth_val <- functional_boxplot(
      boot_c,
      depth_method = method, emp_factor = ef
    )$depth_values
    curve50 <- boot_c[depth_val > sort(depth_val)[dim(boot_c)[1] * 0.5], ]
    curve90 <- boot_c[depth_val > sort(depth_val)[dim(boot_c)[1] * 0.1], ]
    curve95 <- boot_c[depth_val > sort(depth_val)[dim(boot_c)[1] * 0.05], ]
    bounds_c50 <- test_c[m] + cbind(
      apply(curve50, 2, min), apply(curve50, 2, max)
    )
    bounds_c90 <- test_c[m] + cbind(
      apply(curve90, 2, min), apply(curve90, 2, max)
    )
    bounds_c95 <- test_c[m] + cbind(
      apply(curve95, 2, min), apply(curve95, 2, max)
    )
  } else {
    boot_i <- A
    boot_c <- boot_i
    bounds_i50 <- quantile(A, probs = c(0.25, 0.75))
    bounds_i90 <- quantile(A, probs = c(0.05, 0.95))
    bounds_i95 <- quantile(A, probs = c(0.025, 0.975))
    bounds_c50 <- test_c[m] + bounds_i50
    bounds_c90 <- test_c[m] + bounds_i90
    bounds_c95 <- test_c[m] + bounds_i95
  }

  list(
    test_i = test * ep, test_c = test_c * ep,
    pred_i = pred * ep, pred_c = pred_c * ep,
    boot_i = boot_i * ep, boot_c = boot_c * ep,
    bounds_i50 = bounds_i50 * ep, bounds_i90 = bounds_i90 * ep,
    bounds_i95 = bounds_i95 * ep, bounds_c50 = bounds_c50 * ep,
    bounds_c90 = bounds_c90 * ep, bounds_c95 = bounds_c95 * ep
  )
}
