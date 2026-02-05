# nolint: start
library(scales)

plot_bag_covariate <- function(data, covariate, factor = 2.5, depth_method = "extremal") {
  # Initialize list to store hulls and medians
  hulls <- list()
  medians <- list()
  covariate_levels <- levels(data[[covariate]])
  
  # Process each category
  for (i in seq_along(covariate_levels)) {
    temp <- data %>% 
      filter(.data[[covariate]] == covariate_levels[i]) %>% 
      dplyr::select(starts_with("lr_incpaid")) %>% 
      data.matrix()
    
    depth_model <- functional_boxplot(
      temp, depth_method, emp_factor = factor
    )
    
    inner50 <- calculate_depth_groups(depth_model$depth_values, nrow(temp), 0.5)
    
    # Store median curve
    medians[[i]] <- temp[depth_model$median_curve, ] %>%
      as.data.frame() %>%
      pivot_longer(
        cols = everything(),
        names_to = "lag",
        values_to = "value"
      ) %>%
      mutate(
        lag = rep(0:9, length.out = n()),
        group = covariate_levels[i]
      )
    
    # Create hull for this category
    hulls[[i]] <- temp[inner50, ] %>%
      as.data.frame() %>%
      pivot_longer(
        cols = everything(),
        names_to = "lag",
        values_to = "value"
      ) %>%
      mutate(
        lag = rep(0:9, length.out = n()),
        group = covariate_levels[i]
      ) %>%
      group_by(lag) %>%
      summarise(
        min_val = min(value),
        max_val = max(value),
        group = first(group),
        .groups = "drop"
      )
  }
  
  # Combine hulls and medians for plotting
  plot_data <- bind_rows(hulls)
  median_data <- bind_rows(medians)
  
  # Create plot
  plot <- ggplot() +
    # Add upper and lower bands
    geom_line(
      data = plot_data,
      aes(x = lag, y = min_val, color = group),
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_line(
      data = plot_data,
      aes(x = lag, y = max_val, color = group),
      linetype = "dashed",
      linewidth = 1
    ) +
    # Add median curves
    geom_line(
      data = median_data,
      aes(x = lag, y = value, color = group),
      linewidth = 1.5
    ) +
    scale_color_manual(
      name = gsub("_", " ", tools::toTitleCase(covariate)),
      values = setNames(c(
          "purple", "blue", "green", "orange", "red", "brown"
        )[seq_along(covariate_levels)], covariate_levels
      )
    ) +
    labs(
      x = "Lag",
      y = "ILR",
      # title = paste(
      #   "Median Curve and Inner 50% Bands by", 
      #   tools::toTitleCase(gsub("_", " ", covariate))
      # )
    ) +
    theme_bw() +
    theme(
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = "white"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text = element_text(size = 24),
      axis.title = element_text(size = 24),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    scale_x_continuous(breaks = seq(0, 9, 1))
  
  return(plot)
}

plot_median_accident_years <- function(data, factor = 2.5, depth_method = "extremal") {
  highlight_years <- c(1990, 1995, 2000, 2005, 2010)
  medians <- list()
  accident_years <- sort(unique(data$accident_year))
  selected_years <- intersect(accident_years, highlight_years)
  
  for (i in seq_along(selected_years)) {
    temp <- data %>%
      filter(accident_year == selected_years[i]) %>%
      dplyr::select(starts_with("lr_incpaid")) %>%
      data.matrix()
    depth_model <- functional_boxplot(temp, depth_method, emp_factor = factor)
    medians[[i]] <- temp[depth_model$median_curve, ] %>%
      as.data.frame() %>%
      pivot_longer(
        cols = everything(),
        names_to = "lag",
        values_to = "value"
      ) %>%
      mutate(
        lag = rep(0:9, length.out = n()),
        group = as.character(selected_years[i])
      )
  }
  
  median_data <- bind_rows(medians)
  year_colors <- setNames(
    #c("purple", "blue", "green", "orange", "red") # "#E69F00", "#56B4E9", 
    c("#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[seq_along(selected_years)],
    as.character(selected_years)
  )
  
  plot <- ggplot(median_data, aes(x = lag, y = value, color = group)) +
    geom_line(linewidth = 1.75) +
    scale_color_manual(
      name = "Accident Year",
      values = year_colors
    ) +
    labs(
      x = "Lag", y = "ILR",
      # title = "Median ILR Curves for Selected Accident Years"
    ) +
    theme_bw() +
    theme(
      legend.position = c(0.9, 0.95),
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = "white"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text = element_text(size = 24),
      axis.title = element_text(size = 24),
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 24),
      panel.border = element_rect(color = "white", fill = NA, linewidth = 1),
      axis.line.y.left = element_line(color = "black"), # Keep left axis line
      axis.line.x.bottom = element_line(color = "black"), # Keep bottom axis line
      axis.line.y.right = element_blank(), # Remove right axis line
      axis.line.x.top = element_blank() # Remove top axis line
    ) +
    scale_x_continuous(breaks = seq(0, 9, 1)) + 
    coord_cartesian(ylim=c(0,0.28), xlim=c(-0.02,9.02),expand=FALSE)
  
  return(plot)
}

plot_ilr_by_accident_year_and_lag <- function(data) {
  ilr_long <- data %>%
    pivot_longer(
      cols = starts_with("lr_incpaid"),
      names_to = "lag",
      values_to = "ilr"
    ) %>%
    mutate(
      lag = as.numeric(gsub("lr_incpaid\\.", "", lag)),
      accident_year = as.numeric(accident_year)
    )
  
  summary_ilr <- ilr_long %>%
    group_by(accident_year, lag) %>%
    summarise(
      median_ilr = median(ilr, na.rm = TRUE),
      q25 = quantile(ilr, 0.25, na.rm = TRUE),
      q75 = quantile(ilr, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Compute y-limits for each left panel (even lags)
  left_lags <- seq(0, 8, by = 2)
  y_limits <- lapply(left_lags, function(lg) {
    dat <- summary_ilr %>% filter(lag >= lg, lag <= lg + 1)
    c(min(dat$q25, na.rm = TRUE), max(dat$q75, na.rm = TRUE))
  })
  
  # Assign y-limits to both panels in each row
  y_scales <- lapply(1:10, function(i) {
    row_idx <- ((i - 1) %/% 2) + 1
    lims <- y_limits[[row_idx]]
    if (i %% 2 == 1) {
      # Left panel: show axis
      scale_y_continuous(limits = lims, breaks=pretty_breaks(n = 3))
    } else {
      # Right panel: same limits, but hide axis
      scale_y_continuous(limits = lims, labels = NULL, breaks = NULL)
    }
  })
  y_labels = data.frame(lag=0:9, labs=sapply(0:9, function(i) {paste("s=",i,sep="")}), ylims=sapply(1:10, function(i) {
    row_idx <- ((i - 1) %/% 2) + 1
    y_limits[[row_idx]][2]*0.95 }) )
  
  plot <- ggplot(summary_ilr, aes(x = accident_year, y = median_ilr)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "gray80", alpha = 0.5) +
    geom_line(color = "blue", linewidth = 1.5) +
    #annotate("text", x= 2004, y = ylims, label =paste("s="+lag))+
    geom_text(aes(x=1994, y=ylims, label =  labs ), data = y_labels, vjust = 1, size=8) +
    ggh4x::facet_wrap2(~ lag, ncol = 2, scales = "free_y", axes="all", remove_labels = "x") +
    labs(
      x = "Accident Year", y = "ILR",
      # title = "Median and Interquartile Range of ILR by Accident Year and Lag"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text = element_text(size = 24),
      axis.title = element_text(size = 24),
      panel.grid = element_blank(),
      strip.text.x = element_blank(),
      strip.background = element_blank(),
      #strip.text = element_text(size = 16),
      panel.border = element_rect(color = "white", fill = NA, linewidth = 1),
      axis.line.y.left = element_line(color = "black"), # Keep left axis line
      axis.line.x.bottom = element_line(color = "black"), # Keep bottom axis line
      axis.line.y.right = element_blank(), # Remove right axis line
      axis.line.x.top = element_blank() # Remove top axis line
    ) + ggh4x::facetted_pos_scales(y = y_scales) +
    scale_x_continuous(breaks = seq(1990, 2005, 5)) + 
    coord_cartesian(xlim=c(1987,2009),expand=FALSE)
  
  return(plot)
}
