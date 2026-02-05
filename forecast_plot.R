predict_plot <- function(m, company) {
  # Use parameterized depth_interval with global variables
  res <- depth_interval(
    m, company, test_cov_all, test_matrix_all, pred_all, boot_pred_b_all
  )
  bounds_c <- data.frame(
    lag = m:9,
    pred = res$pred_c,
    actual = res$test_c[(m + 1):10],
    lower50 = res$bounds_c50[, 1],
    upper50 = res$bounds_c50[, 2],
    lower90 = res$bounds_c90[, 1],
    upper90 = res$bounds_c90[, 2],
    lower95 = res$bounds_c95[, 1],
    upper95 = res$bounds_c95[, 2]
  )
  bounds_i <- data.frame(
    lag = m:9,
    pred = res$pred_i,
    actual = res$test_i[(m + 1):10],
    lower50 = res$bounds_i50[, 1],
    upper50 = res$bounds_i50[, 2],
    lower90 = res$bounds_i90[, 1],
    upper90 = res$bounds_i90[, 2],
    lower95 = res$bounds_i95[, 1],
    upper95 = res$bounds_i95[, 2]
  )

  ub_i <- max(bounds_i$upper95, bounds_i$actual)
  lb_i <- min(bounds_i$lower95, bounds_i$actual)
  ub_c <- max(bounds_c$upper95, bounds_c$actual)
  lb_c <- min(bounds_c$lower95, bounds_c$actual)

  if (m < 9) {
    plt_i <- ggplot(bounds_i, aes(lag)) +
      geom_ribbon(aes(ymin = lower95, ymax = upper95), fill = "#F0F0F0") +
      geom_ribbon(aes(ymin = lower90, ymax = upper90), fill = "#D4D4D4") +
      geom_ribbon(aes(ymin = lower50, ymax = upper50), fill = "#B3B3B3") +
      geom_line(aes(y = pred), color = "red", linetype = 2, linewidth = 2) +
      geom_line(aes(y = actual), color = "black", linewidth = 2) +
      ylab("Incremental Loss Ratio") +
      theme_bw() +
      theme(
        legend.position = "inside",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      ) +
      scale_x_continuous(limits = c(3, 9), breaks = 3:9) +
      scale_y_continuous(
        name   = "Incremental Losses (million)",
        limits = c(lb_i, ub_i),
        labels = function(x) formatC(x / 1e6, format = "f", digits = 0)
      )

    plt_c <- ggplot(bounds_c, aes(lag)) +
      geom_ribbon(aes(ymin = lower95, ymax = upper95), fill = "#F0F0F0") +
      geom_ribbon(aes(ymin = lower90, ymax = upper90), fill = "#D4D4D4") +
      geom_ribbon(aes(ymin = lower50, ymax = upper50), fill = "#B3B3B3") +
      geom_line(aes(y = pred), color = "red", linewidth = 2, linetype = 2) +
      geom_line(aes(y = actual), color = "black", linewidth = 2) +
      ylab("Cumulative Loss Ratio") +
      theme_bw() +
      theme(
        legend.position = "inside",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      ) +
      scale_x_continuous(limits = c(3, 9), breaks = 3:9) +
      scale_y_continuous(
        name   = "Cumulative Losses (million)",
        limits = c(lb_c, ub_c),
        labels = function(x) formatC(x / 1e6, format = "f", digits = 0)
      )
  } else {
    plt_i <- ggplot(bounds_i, aes(x = lag)) +
      geom_point(aes(y = lower95), color = "#F0F0F0") +
      geom_point(aes(y = upper95), color = "#F0F0F0") +
      geom_point(aes(y = lower90), color = "#D4D4D4") +
      geom_point(aes(y = upper90), color = "#D4D4D4") +
      geom_point(aes(y = lower50), color = "#B3B3B3") +
      geom_point(aes(y = upper50), color = "#B3B3B3") +
      geom_point(aes(y = pred), color = "red", shape = 2) +
      geom_point(aes(y = actual), color = "black") +
      theme_bw() +
      theme(
        legend.position = "inside",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      ) +
      scale_x_continuous(limits = c(3, 9), breaks = 3:9) +
      scale_y_continuous(
        name   = "Incremental Losses (million)",
        limits = c(lb_i, ub_i),
        labels = function(x) formatC(x / 1e6, format = "f", digits = 0)
      )

    plt_c <- ggplot(bounds_c, aes(x = lag)) +
      geom_point(aes(y = lower95), color = "#F0F0F0") +
      geom_point(aes(y = upper95), color = "#F0F0F0") +
      geom_point(aes(y = lower90), color = "#D4D4D4") +
      geom_point(aes(y = upper90), color = "#D4D4D4") +
      geom_point(aes(y = lower50), color = "#B3B3B3") +
      geom_point(aes(y = upper50), color = "#B3B3B3") +
      geom_point(aes(y = pred), color = "red", shape = 2) +
      geom_point(aes(y = actual), color = "black") +
      theme_bw() +
      theme(
        legend.position = "inside",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      ) +
      scale_x_continuous(limits = c(3, 9), breaks = 3:9) +
      scale_y_continuous(
        name   = "Cumulative Losses (million)",
        limits = c(lb_c, ub_c),
        labels = function(x) formatC(x / 1e6, format = "f", digits = 0)
      )
  }

  list(
    lb_i = lb_i, ub_i = ub_i,
    lb_c = lb_c, ub_c = ub_c,
    plt_i = plt_i, plt_c = plt_c
  )
}

predict_cl_plot <- function(company) {
  out <- list()
  for (m in 3:6) {
    # Use parameterized interval functions with global variables
    obj1 <- pls_interval(
      m, company, test_cov_all, test_matrix_all, pred_all, boot_pred_b_all
    )
    obj2 <- depth_interval(
      m, company, test_cov_all, test_matrix_all, pred_all, boot_pred_b_all
    )
    obj3 <- cl(company)
    out[[m - 2]] <- data.frame(
      Lag = m:9, actual = obj1$test_c[(m + 1):10],
      PRED = obj1$pred_c,
      ED.LB = obj2$bounds_c95[, 1],
      ED.UB = obj2$bounds_c95[, 2],
      CLPRED = obj3$pred[[m]],
      CL.LB = obj3$pred[[m]] - qnorm(0.975) * sqrt(obj3$msep[[m]]),
      CL.UB = obj3$pred[[m]] + qnorm(0.975) * sqrt(obj3$msep[[m]])
    )
  }
  ub <- max(out[[1]]$ED.UB, out[[1]]$CL.UB, out[[1]]$actual)
  lb <- min(out[[1]]$ED.LB, out[[1]]$CL.LB, out[[1]]$actual)
  a <- (ub - lb) / 2

  plt <- list()
  for (i in 1:4) {
    if (i > 1) {
      ub <- out[[i]]$PRED[8 - i] + a
      lb <- out[[i]]$PRED[8 - i] - a
    }
    plt[[i]] <- ggplot(out[[i]], aes(x = Lag)) +
      # Gray shaded region
      geom_ribbon(aes(ymin = ED.LB, ymax = ED.UB, fill = "EXD Central Region"), alpha = 0.5) +
      # Lines with proper aes mapping
      geom_line(aes(y = PRED, color = "PLS Prediction"), linetype = 2, linewidth = 1) +
      geom_line(aes(y = actual, color = "Actual"), linewidth = 1) +
      geom_line(aes(y = CLPRED, color = "CL Prediction and Interval"), linetype = 4, linewidth = 1) +
      geom_line(aes(y = CL.LB, color = "CL Prediction and Interval"), linetype = 3, linewidth = 1) +
      geom_line(aes(y = CL.UB, color = "CL Prediction and Interval"), linetype = 3, linewidth = 1) +
      # Custom colors and fill for legend
      scale_color_manual(
        name = "",
        values = c(
          "PLS Prediction" = "red",
          "Actual" = "black",
          "CL Prediction and Interval" = "blue"
        )
      ) +
      scale_fill_manual(name = "", values = c("EXD Central Region" = "#D4D4D4")) +
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
      scale_x_continuous(limits = c(3, 9), breaks = 3:9) +
      scale_y_continuous(
        name   = "Cumulative Losses (million)",
        limits = c(lb, ub),
        labels = function(x) formatC(x / 1e6, format = "f", digits = 0)
      )
    # Order legends explicitly: color first, then fill
    # + guides(color = guide_legend(order = 1),
    #        fill = guide_legend(order = 2))
  }
  plt
}

# =============================================================================
# FIXED ORIGIN BACKTEST PLOTTING FUNCTIONS
# =============================================================================
# These functions are used for Figures 5 and 9 in both main.R and replication.R
# They require global variables: test_cov_all_fix, test_matrix_all_fix,
# pred_all_fix, boot_pred_b_all_fix
# Note: depth_interval_fix() is in prediction_interval.R

library(patchwork)

#' Prediction plot for fixed origin backtest
#' @param m Number of completed lags
#' @param company Company SNL key
#' @param method Depth method: "extremal", "mbd", "tvd", "bd"
#' @return List with incremental plot (plt_i), cumulative plot (plt_c), and data
predict_plot_fix <- function(m, company, method = "bd") {
  res <- depth_interval_fix(m, company, method = method)

  if (m < 9) {
    prediction_m <- data.frame(
      ILR = as.numeric(
        c(
          res$test_i, res$test_i[m], res$pred_i, res$test_i[m],
          res$bounds_i95[, 1], res$test_i[m], res$bounds_i95[, 2]
        )
      ),
      Lag = c(0:9, (m - 1):9, (m - 1):9, (m - 1):9),
      Group = c(rep(1, 10), rep(2, 11 - m), rep(3, 11 - m), rep(4, 11 - m))
    )
    prediction_mc <- data.frame(
      CLR = as.numeric(
        c(
          res$test_c, res$test_c[m], res$pred_c, res$test_c[m],
          res$bounds_c95[, 1], res$test_c[m], res$bounds_c95[, 2]
        )
      ),
      Lag = c(0:9, (m - 1):9, (m - 1):9, (m - 1):9),
      Group = c(rep(1, 10), rep(2, 11 - m), rep(3, 11 - m), rep(4, 11 - m))
    )
  } else {
    prediction_m <- data.frame(
      ILR = as.numeric(
        c(
          res$test_i, res$test_i[m], res$pred_i, res$test_i[m],
          res$bounds_i95[1], res$test_i[m], res$bounds_i95[2]
        )
      ),
      Lag = c(0:9, (m - 1):9, (m - 1):9, (m - 1):9),
      Group = c(rep(1, 10), rep(2, 11 - m), rep(3, 11 - m), rep(4, 11 - m))
    )
    prediction_mc <- data.frame(
      CLR = as.numeric(
        c(
          res$test_c, res$test_c[m], res$pred_c, res$test_c[m],
          res$bounds_c95[1], res$test_c[m], res$bounds_c95[2]
        )
      ),
      Lag = c(0:9, (m - 1):9, (m - 1):9, (m - 1):9),
      Group = c(rep(1, 10), rep(2, 11 - m), rep(3, 11 - m), rep(4, 11 - m))
    )
  }

  plt_i <- ggplot(prediction_m, aes(x = Lag, y = ILR, linetype = factor(Group))) +
    geom_line(size = 1.5) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 24),
      axis.title = element_text(size = 24),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    scale_x_continuous(breaks = 0:9)

  plt_c <- ggplot(prediction_mc, aes(x = Lag, y = CLR, linetype = factor(Group))) +
    geom_line(size = 1.5) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 24),
      axis.title = element_text(size = 24),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    scale_x_continuous(breaks = 0:9)

  list(plt_i = plt_i, plt_c = plt_c, df = prediction_mc)
}

#' Create grid of prediction plots for all lags (m = 1 to 9)
#' @param company Company SNL key
#' @param type "i" for incremental, "c" for cumulative
#' @param lb Lower bound for y-axis
#' @param ub Upper bound for y-axis
#' @param method Depth method
#' @return Combined plot grid
predict_plot_grid_fix <- function(company, type = "c", lb = NULL, ub = NULL,
                                  method = "bd") {
  plot_list <- lapply(1:9, function(i) {
    pp <- predict_plot_fix(m = i, company = company, method = method)
    p <- if (type == "i") pp$plt_i else pp$plt_c

    if (!is.null(lb) && !is.null(ub)) {
      p <- p + ylim(lb, ub)
    }
    p <- p + theme(plot.margin = margin(2, 2, 2, 2))

    # Adjust axis labels based on position in 3x3 grid
    if (i %in% c(1, 4)) {
      p + theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    } else if (i %in% c(8, 9)) {
      p + theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    } else if (i == 7) {
      p
    } else {
      p + theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )
    }
  })

  wrap_plots(plot_list, ncol = 3)
}

#' Plot CLR prediction and interval tracking as more lags become known
#' @param company Company SNL key
#' @param method Depth method
#' @return ggplot object
predict_clr_tracking_fix <- function(company, method = "bd") {
  prediction_clr <- NULL

  for (m in 1:9) {
    for (g in 1:4) {
      pp <- predict_plot_fix(m, company, method = method)$df %>%
        filter(Group == g)
      if (g == 1) {
        prediction_clr <- prediction_clr %>%
          rbind(data.frame(CLR = pp[10, 1], Group = g, m = m))
      } else {
        prediction_clr <- prediction_clr %>%
          rbind(data.frame(CLR = pp[11 - m, 1], Group = g, m = m))
      }
    }
  }

  ggplot(
    prediction_clr,
    aes(x = m, y = CLR, color = factor(Group), linetype = factor(Group))
  ) +
    geom_line(size = 1.5) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 24),
      axis.title = element_text(size = 24),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    scale_x_continuous(breaks = 0:9) +
    xlab("Lags Completed") +
    ylab("Ultimate Cumulative Losses") +
    scale_color_manual(values = c("black", "red", "blue", "blue")) +
    scale_linetype_manual(values = c(1, 2, 4, 4))
}
