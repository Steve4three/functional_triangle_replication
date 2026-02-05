# nolint: start
# Read triangle and company info data and combine.
# Also product data matrix for PCA.
read_data <- function(triangle_path, company_path) {
  data <- read.csv(triangle_path)
  time_invariant_vars <- read.csv(company_path)
  colnames(time_invariant_vars)[1] <- "snl_key"
  data %>% full_join(time_invariant_vars, by = "snl_key")
}

process_loss_ratios <- function(data) {
  data %>%
    filter(accident_year == data_year - 9) %>%
    mutate(
      lr_paid = paid_lr / 100,
      lr_reported = reported_lr / 100,
      lr_incurred = incurred_lr / 100,
      ep = ep * 1000
    ) %>%
    mutate(
      loss_paid = lr_paid * ep,
      loss_reported = lr_reported * ep,
      loss_incurred = lr_incurred * ep
    ) %>%
    dplyr::select(-paid_lr, -reported_lr, -incurred_lr)
}

reshape_and_calculate_ldf <- function(data) {
  reshape(data,
    v.names = c("loss_paid", "loss_reported", "loss_incurred",
                "lr_paid", "lr_reported", "lr_incurred"),
    idvar = c("snl_key", "lob", "accident_year", "ep"),
    timevar = "dev_lag",
    direction = "wide"
  ) %>%
    mutate(
      ldf_paid.1 = loss_paid.1 / loss_paid.0,
      ldf_paid.2 = loss_paid.2 / loss_paid.1,
      ldf_paid.3 = loss_paid.3 / loss_paid.2,
      ldf_paid.4 = loss_paid.4 / loss_paid.3,
      ldf_paid.5 = loss_paid.5 / loss_paid.4,
      ldf_paid.6 = loss_paid.6 / loss_paid.5,
      ldf_paid.7 = loss_paid.7 / loss_paid.6,
      ldf_paid.8 = loss_paid.8 / loss_paid.7,
      ldf_paid.9 = loss_paid.9 / loss_paid.8,
      lr_incpaid.0 = lr_paid.0,
      lr_incpaid.1 = lr_paid.1 - lr_paid.0,
      lr_incpaid.2 = lr_paid.2 - lr_paid.1,
      lr_incpaid.3 = lr_paid.3 - lr_paid.2,
      lr_incpaid.4 = lr_paid.4 - lr_paid.3,
      lr_incpaid.5 = lr_paid.5 - lr_paid.4,
      lr_incpaid.6 = lr_paid.6 - lr_paid.5,
      lr_incpaid.7 = lr_paid.7 - lr_paid.6,
      lr_incpaid.8 = lr_paid.8 - lr_paid.7,
      lr_incpaid.9 = lr_paid.9 - lr_paid.8
    )
}

standardize_business_focus <- function(data) {
  data %>%
    mutate(
      business_focus = case_when(
        business_focus == "Commercial General Liability Focus" ~ "Commercial",
        business_focus == "Commercial Lines Focus" ~ "Commercial",
        business_focus == "Commercial Property Focus" ~ "Commercial",
        business_focus == "Commercial Workers Compensation Focus" ~ "WkComp",
        business_focus == "Personal Lines Focus" ~ "Personal",
        TRUE ~ business_focus
      )
    ) %>%
    mutate(business_focus = factor(business_focus,
      levels = c("Commercial", "Personal", "WkComp")
    ))
}

standardize_ownership <- function(data) {
  data %>%
    mutate(
      naic_ownership_structure = case_when(
        naic_ownership_structure == "Mutual Company" ~ "Mutual",
        naic_ownership_structure == "Stock Company" ~ "Stock",
        naic_ownership_structure == "Reciprocal Exchange" ~ "Other",
        TRUE ~ naic_ownership_structure
      )
    ) %>%
    mutate(naic_ownership_structure = factor(naic_ownership_structure,
      levels = c("Mutual", "Stock", "Other")
    ))
}

standardize_geographic_focus <- function(data) {
  data %>%
    mutate(
      geographic_focus = case_when(
        geographic_focus == "Regional - Midwestern Quadrant" ~ "Midwest",
        geographic_focus == "Regional - Northeastern Quadrant" ~ "Northeast",
        geographic_focus == "Regional - Southern Quadrant" ~ "South",
        geographic_focus == "Regional - Western Quadrant" ~ "West",
        TRUE ~ geographic_focus
      )
    ) %>%
    mutate(geographic_focus = factor(geographic_focus,
      levels = c("National", "Midwest", "Northeast", "South", "West")
    ))
}

standardize_fhlb_member <- function(data) {
  data %>%
    mutate(
      fhlb_member_bank = if_else(
        is.na(fhlb_member_bank) | fhlb_member_bank == "",
        "Unknown",
        fhlb_member_bank
      )
    ) %>%
    mutate(fhlb_member_bank = factor(fhlb_member_bank,
      levels = c(
        "Unknown", 
        "Federal Home Loan Bank of Chicago",
        "Federal Home Loan Bank of Dallas", 
        "Federal Home Loan Bank of Des Moines",
        "Federal Home Loan Bank of Indianapolis", 
        "Federal Home Loan Bank of New York"
      )
    ))
}

prepare_data <- function(data, LOB, type = "incremental") {
  data %>%
    filter(lob == LOB) %>%
    dplyr::select(
      snl_key, ep, business_focus, naic_ownership_structure,
      geographic_focus, accident_year, fhlb_member_bank,
      starts_with(if_else(type == "incremental", "lr_incpaid", "lr_paid"))
    ) %>%
    filter(
      business_focus != "Reinsurance",
      business_focus != "Large Reinsurance Focus",
      !(is.na(business_focus)),
      !(is.na(naic_ownership_structure)),
      !(is.na(geographic_focus))
    ) %>%
    mutate(
      id = as.character(
        interaction(accident_year, snl_key)
      )
    ) %>%
    arrange(id)
}

create_pca_matrix <- function(data, type = "incremental") {
  matrix_data <- data.matrix(
    data %>% 
      dplyr::select(
        starts_with(
          if_else(type == "incremental", "lr_incpaid", "lr_paid")
        )
      )
  )
  rownames(matrix_data) <- data$id
  colnames(matrix_data) <- 0:9
  matrix_data
}
# nolint: end