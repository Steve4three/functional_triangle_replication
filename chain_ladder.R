cl <- function(company) {
  # company <- "C3856"
  # Choose one company and filter the company's ILR triangle and covariates
  data_cl <- data_cov %>% filter(snl_key == company)

  ilr <- data_cl %>% dplyr::select(starts_with("lr_incpaid")) %>% as.matrix()
  clr <- t(apply(ilr, 1, cumsum)) # compute CLR
  inc_loss <- data_cl$ep * ilr # Not used for CL factor, but for other methods
  cum_loss <- data_cl$ep * clr
  n_rows <- nrow(cum_loss)

  # Compute CL factors
  f <- c()
  sig_sq <- c()
  for (m in 1:9) {
    idx <- n_rows - m # row of accident year 2010-m
    # chain ladder factor formula
    f[m] <- sum(cum_loss[1:idx, m + 1]) / sum(cum_loss[1:idx, m])
    if (idx > 1) {
      # chain ladder variance factor formula
      sig_num <- sum(
        cum_loss[1:idx, m] *
          (cum_loss[1:idx, m + 1] / cum_loss[1:idx, m] - f[m])^2
      )
      sig_sq[m] <- sig_num / (idx - 1)
    } else {
      sig_sq[m] <- min(
        sig_sq[m - 1],
        sig_sq[m - 2],
        sig_sq[m - 1]^2 / sig_sq[m - 2]
      )
    }
  }
  # Calculate chain ladder prediction and msep
  pred <- list()
  reserve <- c()
  msep <- list()
  for (m in 1:9) {
    idx <- n_rows - m
    # chain ladder cumulative losses prediction for lag m+1 to 10
    pred[[m]] <- cum_loss[idx + 1, m] * cumprod(f[m:9])
    reserve[m] <- cum_loss[idx + 1, m] * (prod(f[m:9]) - 1)
    # chain ladder mse of prediction of cumulative loss (Theorem 3, Mack 1993)
    c1 <- ifelse(
      m < 9,
      c(cum_loss[idx + 1, m], pred[[m]][1:(9 - m)]),
      cum_loss[idx + 1, m]
    )
    s1 <- c()
    for (i in 1:(10 - m)) {
      s1[i] <- sum(cum_loss[1:(idx + 1 - i), i])
    }
    msep_part1 <- pred[[m]]^2 * cumsum(sig_sq[m:9] / f[m:9]^2 / c1)
    msep_part2 <- pred[[m]]^2 * cumsum(sig_sq[m:9] / f[m:9]^2 / s1)
    msep[[m]] <- msep_part1 + msep_part2
  }
  
  list(pred = pred, reserve = reserve, msep = msep)
}

# Function to calculate reserves for a single company
calculate_reserves <- function(company) {
  # Get chain ladder results
  cl_res_company <- sum(cl(company)$reserve)
  
  # Get data for the company
  data_cl <- data_cov %>% filter(snl_key == company)
  ilr <- data_cl %>% dplyr::select(starts_with("lr_incpaid")) %>% as.matrix()
  clr <- t(apply(ilr, 1, cumsum))
  cum_loss <- data_cl$ep * clr
  n_rows <- nrow(cum_loss)
  
  # Initialize results for this company
  pls_res_company <- 0
  tru_res_company <- 0
  cl_err_company <- rep(0, 9)
  pls_err_company <- rep(0, 9)
  
  # Calculate reserves for each lag
  for (m in 1:9) {
    # Chain ladder reserve
    cl_res_m <- cl(company)$reserve[m]
    
    # PLS reserve
    idx <- which(test_cov_all[[m]]$snl_key == company)
    pls_res_m <- sum(pred_all[[m]][idx, ]) * data_cl$ep[n_rows - m + 1]
    pls_res_company <- pls_res_company + pls_res_m
    
    # True reserve
    tru_res_m <- cum_loss[n_rows - m + 1, 10] - cum_loss[n_rows - m + 1, m]
    tru_res_company <- tru_res_company + tru_res_m
    
    # Errors
    cl_err_company[m] <- abs(cl_res_m - tru_res_m)
    pls_err_company[m] <- abs(pls_res_m - tru_res_m)
  }
  
  list(
    cl_res = cl_res_company,
    pls_res = pls_res_company,
    tru_res = tru_res_company,
    cl_err = cl_err_company,
    pls_err = pls_err_company
  )
} 
