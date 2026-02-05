# nolint: start
# p: percentage of training
# fold: TRUE if k-fold data
# fn: number of folds
# n: fold n as test set
data_gen <- function(data, cutoff1, cutoff2) {
  train_cov <- data %>%
    filter(accident_year <= cutoff1) %>%
    dplyr::select(
      snl_key, ep, accident_year,
      business_focus, naic_ownership_structure,
      geographic_focus
    ) %>%
    mutate(time = accident_year - 1986)

  train_matrix <- data %>%
    filter(accident_year <= cutoff1) %>%
    dplyr::select(starts_with("LR_incpaid")) %>%
    data.matrix()
  colnames(train_matrix) <- 0:9

  test_cov <- data %>%
    filter(
      accident_year > cutoff1,
      accident_year <= cutoff2
    ) %>%
    dplyr::select(
      snl_key, ep, accident_year,
      business_focus, naic_ownership_structure,
      geographic_focus
    ) %>%
    mutate(time = accident_year - 1986)

  test_matrix <- data %>%
    filter(
      accident_year > cutoff1,
      accident_year <= cutoff2
    ) %>%
    dplyr::select(starts_with("LR_incpaid")) %>%
    data.matrix()
  colnames(test_matrix) <- 0:9

  # Run PCA on training data and obtain 10 x 1 median curve 10 x n loading
  # matrix and nobs(train) x n betas

  pca <- prcomp(train_matrix, 10)
  # Create dataframe with betas and covariates for beta forecast
  beta <- pca$x
  colnames(beta) <- paste0("beta", 1:10)

  list(
    train_cov = train_cov,
    test_cov = test_cov,
    train_matrix = train_matrix,
    test_matrix = test_matrix,
    pca = pca,
    beta = beta
  )
}

data_split <- function(data, p = NULL) {
  partition <- createDataPartition(data$depth, p = p, list = FALSE)

  train_cov <- data[partition, ] %>%
    dplyr::select(
      snl_key, ep, accident_year,
      business_focus, naic_ownership_structure,
      geographic_focus
    ) %>%
    mutate(time = accident_year - 1986)

  test_cov <- data[-partition, ] %>%
    dplyr::select(
      snl_key, ep, accident_year,
      business_focus, naic_ownership_structure,
      geographic_focus
    ) %>%
    mutate(time = accident_year - 1986)

  train_matrix <- data[partition, ] %>%
    dplyr::select(starts_with("LR_incpaid")) %>%
    data.matrix()
  colnames(train_matrix) <- 0:9

  test_matrix <- data[-partition, ] %>%
    dplyr::select(starts_with("LR_incpaid")) %>%
    data.matrix()
  colnames(test_matrix) <- 0:9

  pca <- prcomp(train_matrix, 10)
  # Create dataframe with betas and covariates for beta forecast
  beta <- pca$x

  colnames(beta) <- paste0("beta", 1:10)

  list(
    train_cov = train_cov,
    test_cov = test_cov,
    train_matrix = train_matrix,
    test_matrix = test_matrix,
    pca = pca,
    beta = beta
  )
}

fold_gen <- function(data, fn = 5) {
  train_cov <- list()
  test_cov <- list()
  train_matrix <- list()
  test_matrix <- list()
  pca <- list()
  beta <- list()
  folds <- createFolds(data$depth, k = fn, list = TRUE)

  for (n in 1:fn) {
    train_cov[[n]] <- data[-folds[[n]], ] %>%
      dplyr::select(
        snl_key, ep, accident_year,
        business_focus, naic_ownership_structure,
        geographic_focus
      ) %>%
      mutate(time = accident_year - 1986)

    test_cov[[n]] <- data[folds[[n]], ] %>%
      dplyr::select(
        snl_key, ep, accident_year,
        business_focus, naic_ownership_structure,
        geographic_focus
      ) %>%
      mutate(time = accident_year - 1986)

    train_matrix[[n]] <- data[-folds[[n]], ] %>%
      dplyr::select(starts_with("LR_incpaid")) %>%
      data.matrix()
    colnames(train_matrix[[n]]) <- 0:9

    test_matrix[[n]] <- data[folds[[n]], ] %>%
      dplyr::select(starts_with("LR_incpaid")) %>%
      data.matrix()
    colnames(test_matrix[[n]]) <- 0:9

    pca[[n]] <- prcomp(train_matrix[[n]], 10)
    # Create dataframe with betas and covariates for beta forecast
    beta[[n]] <- pca[[n]]$x

    colnames(beta[[n]]) <- paste0("beta", 1:10)
  }

  list(
    train_cov = train_cov,
    test_cov = test_cov,
    train_matrix = train_matrix,
    test_matrix = test_matrix,
    pca = pca,
    beta = beta
  )
}

beta_reg <- function(beta_data, covariate, k, alpha = NULL, lasso = TRUE) {
  if (lasso == FALSE) {
    model <- NULL
    for (n in 1:k) {
      x <- model.matrix(~ ., data = covariate %>%
        mutate(
          log_ep = log(ep),
          log_ep_time = log(ep) * time
        ) %>%
        dplyr::select(-snl_key, -ep, -accident_year)
      ) %>%
        data.frame() %>%
        dplyr::select(-X.Intercept.) %>%
        add_column(beta = beta_data[, n])

      model_full <- glm(
        beta ~ .,
        family = gaussian(),
        weights = log_ep,
        data = x
      )
      model[[n]] <- stepAIC(
        model_full,
        k = log(nrow(beta_data)),
        trace = FALSE
      )
    }
  } else {
    x <- model.matrix(~ ., data = covariate %>%
      mutate(
        log_ep = log(ep),
        log_ep_time = log(ep) * time
      ) %>%
      dplyr::select(-snl_key, -ep, -accident_year)
    )

    model <- NULL
    for (n in 1:k) {
      cv_lam <- cv.glmnet(
        x = x,
        y = beta_data[, n],
        family = "gaussian",
        alpha = alpha
      )
      model[[n]] <- glmnet(
        x = x,
        y = beta_data[, n],
        family = "gaussian",
        weights = x[, "log_ep"],
        alpha = alpha,
        lambda = cv_lam$lambda.min
      )
    }
  }
  model
}

pls_triangle <- function(
  lambda, m, beta_model, covariate, test_matrix, pca, k, lasso = TRUE
) {
  center <- pca$center
  loading <- pca$rotation

  x <- covariate %>%
    mutate(
      log_ep = log(ep),
      log_ep_time = log(ep) * time
    ) %>%
    dplyr::select(-snl_key, -ep, -accident_year) %>%
    model.matrix(~ ., data = .) %>%
    data.frame() 
  
  if (lasso == TRUE) {
    x <- x %>% as.matrix()
  } else {
    x <- x %>% dplyr::select(-X.Intercept.)
  }

  beta_pred <- sapply(1:k, function(n) {
    predict(
      beta_model[[n]],
      x,
      type = "response"
    )
  })

  # Extract yc for all rows at once
  y_c <- test_matrix[, 1:m] - matrix(
    center[1:m], nrow(test_matrix), m, byrow = TRUE
  )

  # Compute beta_pls for all rows
  fe <- matrix(loading[1:m, 1:k], m, k)
  fl <- matrix(loading[(m + 1):10, 1:k], 10 - m, k)
  i_k <- diag(k)
  fe_t_fe <- solve(t(fe) %*% fe + lambda * i_k)
  fe_t_yc <- t(fe) %*% t(y_c) + lambda * t(beta_pred)
  beta_pls <- fe_t_fe %*% fe_t_yc

  # Calculate predictions for all rows
  pred <- matrix(
    center[(m + 1):10], nrow(test_matrix), 10 - m, byrow = TRUE
  ) + t(fl %*% beta_pls)

  list(pred = pred, beta_pred = beta_pred, beta_pls = t(beta_pls))
}

# nolint: end
