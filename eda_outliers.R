#' Calculate summary statistics for development lag data
#'
#' @param data Matrix of development lag data
#' @param cum Logical indicating if data is cumulative
#'   (removes count statistics)
#' @return Data frame with summary statistics for each lag
#' @export
devlag_summary <- function(data, cum = FALSE) {
  # Initialize summary matrix
  summary_stats <- matrix(0, 15, 10)

  # Define statistic names
  stat_names <- c(
    "Mean", "Std. Dev.", "Median", "MAD", "Min", "Max",
    "95% quantile", "99% quantile", "Skewness", "Kurtosis",
    "5% winsorized mean", "5% winsorized sd",
    "# Positive", "# Zero", "# Negative"
  )

  # Calculate statistics for each lag
  for (i in seq_len(10)) {
    summary_stats[1, i] <- mean(data[, i], na.rm = TRUE) %>% round(4)
    summary_stats[2, i] <- sd(data[, i], na.rm = TRUE) %>% round(4)
    summary_stats[3, i] <- median(data[, i], na.rm = TRUE) %>% round(4)
    summary_stats[4, i] <- mad(data[, i], na.rm = TRUE) %>% round(4)
    summary_stats[5, i] <- min(data[, i], na.rm = TRUE) %>% round(4)
    summary_stats[6, i] <- max(data[, i], na.rm = TRUE) %>% round(4)
    summary_stats[7, i] <- quantile(data[, i], 0.05, na.rm = TRUE) %>% round(4)
    summary_stats[8, i] <- quantile(data[, i], 0.95, na.rm = TRUE) %>% round(4)
    summary_stats[9, i] <- skewness(data[, i], na.rm = TRUE) %>% round(4)
    summary_stats[10, i] <- kurtosis(data[, i], na.rm = TRUE) %>% round(4)
    summary_stats[11, i] <- winsor.mean(data[, i], 0.05) %>% round(4)
    summary_stats[12, i] <- winsor.sd(data[, i], 0.05) %>% round(4)
    summary_stats[13, i] <- sum(data[, i] > 0, na.rm = TRUE) %>% round(0)
    summary_stats[14, i] <- sum(data[, i] == 0, na.rm = TRUE) %>% round(0)
    summary_stats[15, i] <- sum(data[, i] < 0, na.rm = TRUE) %>% round(0)
  }

  # Create result data frame
  result <- cbind(
    data.frame(Statistic = stat_names),
    data.frame(summary_stats)
  )
  colnames(result) <- c("Statistic", 0:9)

  # Remove count statistics for cumulative data
  if (cum) {
    result <- result[-c(13:15), ]
  }

  result
}

#' Calculate depth threshold and inner group
#' @param depth_values Vector of depth values
#' @param n_samples Number of samples
#' @param threshold_pct Percentage threshold (0-1) for inner group
#' @return List containing threshold and inner group indices
calculate_depth_groups <- function(depth_values, n_samples, threshold_pct = 0.5) {
  ind <- data.frame(
    id = seq_len(n_samples),
    dep = depth_values
  )
  threshold <- sort(ind$dep, decreasing = TRUE)[n_samples * threshold_pct]
  
  ind %>% filter(dep >= threshold) %>% pull(id)
}

plot_bag_outliers <- function(pca_data, pca_scores, factor = 2.5, depth_method = "extremal",
                              show.leg="none") {
  # Calculate depth and identify groups
  # show.leg controls the appearance/location of a legend for fill
  depth_model <- functional_boxplot(
    pca_scores, depth_method, emp_factor = factor
  )

  inner50 <- calculate_depth_groups(depth_model$depth_values, nrow(pca_data), 0.5)
  inner90 <- calculate_depth_groups(depth_model$depth_values, nrow(pca_data), 0.9)
  inner95 <- calculate_depth_groups(depth_model$depth_values, nrow(pca_data), 0.95)
  
  # Prepare data for plotting
  plot_data <- pca_data %>%
    as.data.frame() %>%
    mutate(
      id = row_number(),
      group = case_when(
        id %in% depth_model$outliers ~ "Outlier",
        id %in% inner50 ~ "Inner 50%",
        id %in% inner90 ~ "Inner 90%",
        id %in% inner95 ~ "Inner 95%",
        TRUE ~ "Other"
      )
    ) %>%
    pivot_longer(
      cols = -c(id, group),
      names_to = "lag",
      values_to = "value"
    ) %>%
    mutate(lag = as.numeric(lag))
  
  # Create convex hull for inner50, inner90, inner95
  inner50_hull <- plot_data %>%
    filter(group == "Inner 50%") %>%
    group_by(lag) %>%
    summarise(
      min_val = min(value),
      max_val = max(value)
    )
  inner90_hull <- plot_data %>%
    filter(group == "Inner 90%") %>%
    group_by(lag) %>%
    summarise(
      min_val = min(value),
      max_val = max(value)
    )
  inner95_hull <- plot_data %>%
    filter(group == "Inner 95%") %>%
    group_by(lag) %>%
    summarise(
      min_val = min(value),
      max_val = max(value)
    )
  inner50_hull$maxs <- sapply(1:10, function(i) {max(pca_data_i[inner50,i])})
  inner90_hull$maxs <- sapply(1:10, function(i) {max(pca_data_i[inner90,i])})
  inner95_hull$maxs <- sapply(1:10, function(i) {max(pca_data_i[inner95,i])})
  
  inner50_hull$mins <- sapply(1:10, function(i) {min(pca_data_i[inner50,i])})
  inner90_hull$mins <- sapply(1:10, function(i) {min(pca_data_i[inner90,i])})
  inner95_hull$mins <- sapply(1:10, function(i) {min(pca_data_i[inner95,i])})
  median_curve <- plot_data %>% filter(id == depth_model$median_curve)
  # Create plot
  plot <- ggplot(plot_data, aes(x = lag, y = value, group = id)) +
    # Add outlier curves as dotted lines
    geom_line(
      data = filter(plot_data, group == "Outlier"),
      aes(linetype = "Outlier", color = "Outlier"),
      color = "black",
      linetype = "dotted",
      linewidth = 0.5,
      alpha = 0.5
    ) +
    # Add convex hull for inner95 (lightest)
    geom_ribbon(
      data = inner95_hull,
      aes(x = lag, ymin = mins, ymax = maxs, fill = "Inner 95%"),
      #fill = "#F0F0F0",
      alpha = 1,
      inherit.aes = FALSE
    ) +
    # Add convex hull for inner90 (medium)
    geom_ribbon(
      data = inner90_hull,
      aes(x = lag, ymin = mins, ymax = maxs, fill = "Inner 90%"),
      #fill = "#D4D4D4",
      alpha = 0.5,
      inherit.aes = FALSE
    ) +
    # Add convex hull for inner50 (darkest)
    geom_ribbon(
      data = inner50_hull,
      aes(x = lag, ymin = mins, ymax = maxs, fill = "Inner 50%"),
      #fill = "#B3B3B3",
      alpha = 0.5,
      inherit.aes = FALSE
    ) +
    # Add median curve
    geom_line(
      data = median_curve,
      aes(color = "Median"),
      color = "black",
      linewidth = 1.25
    ) +
    scale_fill_manual(
      name = "Functional Bag",
      values = c(
        "Inner 95%" = "#F0F0F0",
        "Inner 90%" = "#C0C0C0",
        "Inner 50%" = "#A2A2A2"
      ),
      guide = guide_legend()
    ) +
    scale_color_manual(
      name = "Curve",
      values = c(
        "Outlier" = "black",
        "Median" = "black"
      ),
      guide = guide_legend(order = 2)
    ) +
    scale_linetype_manual(
      name = "Curve",
      values = c(
        "Outlier" = "dotted",
        "Median" = "solid"
      ),
      guide = guide_legend(order = 2)
    ) +
    labs(
      x = "Lag", y = "ILR"
    ) +
    theme_bw() +
    theme(
      legend.position = show.leg,
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text = element_text(size = 24),
      axis.title = element_text(size = 24),
      axis.line.y.left = element_line(color = "black"), # Keep left axis line
      axis.line.x.bottom = element_line(color = "black"), # Keep bottom axis line
      axis.line.y.right = element_blank(), # Remove right axis line
      axis.line.x.top = element_blank(), # Remove top axis line
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      panel.border = element_rect(color = "white", fill = NA, linewidth = 1)
    ) +
    scale_x_continuous(breaks = seq(0, 9, 1)) +
    coord_cartesian(xlim=c(0.4, 8.8))
  return(list(
    n = length(depth_model$outliers),
    data = pca_data[-depth_model$outliers, ],
    in50 = pca_data[inner50, ],
    mc = median_curve,
    plt = plot
  ))
}