# =============================================================================
# REPLICATION SCRIPT
# =============================================================================
# This script reproduces all tables and figures from the paper using
# pre-computed results. No raw company data is required.
#
# Required files:
#   - replication_data.Rdata    (3 companies for plotting: P53, C260, P1406)
#   - replication_results.Rdata (aggregate results for tables)
#   - sensitivity_K_results.Rdata (sensitivity analysis results)
#
# Also requires the function files:
#   - functions.R
#   - tuning.R
#   - chain_ladder.R
#   - prediction_interval.R
#   - evaluation.R
#   - forecast_plot.R
# =============================================================================

library(tidyverse)
library(gridExtra)
library(fdaoutlier)

# Source function files (methodology code, no data)
source("./functions.R")
source("./tuning.R")
source("./chain_ladder.R")
source("./prediction_interval.R")
source("./evaluation.R")
source("./forecast_plot.R")

# =============================================================================
# LOAD REPLICATION DATA
# =============================================================================

cat("Loading replication data files...\n")
load("replication_data.Rdata") # 3 companies for plots
load("replication_results.Rdata") # Aggregate results for tables
load("sensitivity_K_results.Rdata") # Sensitivity analysis results
cat("Data loaded successfully.\n")

# =============================================================================
# TABLE 7: MAPE Comparison - All K Values (from sensitivity_K.R)
# =============================================================================

cat("\n=============================================================================\n")
cat("TABLE 7: MAPE Comparison - All K Values\n")
cat("=============================================================================\n")

# Build MAPE comparison table
# Columns: s, K=1, K=2, ..., K=9

mape_table <- data.frame(s = 1:9)

# Add columns for each K value (K=1 to K=9) from sensitivity analysis
for (k in 1:9) {
    k_char <- as.character(k)
    if (k_char %in% names(all_results_sens)) {
        mape_vals <- all_results_sens[[k_char]]$mape
        mape_table[[paste0("K=", k)]] <- ifelse(
            is.na(mape_vals), "-",
            as.character(round(mape_vals, 4))
        )
    } else {
        mape_table[[paste0("K=", k)]] <- "-"
    }
}

print(mape_table)

# Helper function to safely extract coverage/score values
safe_extract <- function(coverage_obj, field_name, col_name) {
    result <- rep(NA, 9)
    if (!is.null(coverage_obj) && !is.null(coverage_obj[[field_name]])) {
        vals <- coverage_obj[[field_name]][, col_name]
        if (length(vals) == 9) {
            result <- vals
        }
    }
    result
}

# =============================================================================
# TABLE 8: Ultimate Coverage Table (EXD)
# =============================================================================

ultimate_coverage_table <- data.frame(s = 1:9)
for (k in 1:9) {
    k_char <- as.character(k)
    if (k_char %in% names(all_results_sens)) {
        cov_vals <- safe_extract(all_results_sens[[k_char]]$coverage, "ultimate_coverage_pct", "EXD")
        ultimate_coverage_table[[paste0("K=", k)]] <- ifelse(
            is.na(cov_vals), "-",
            as.character(round(cov_vals, 4))
        )
    } else {
        ultimate_coverage_table[[paste0("K=", k)]] <- "-"
    }
}
cat("\nUltimate Coverage (EXD) by K:\n")
print(ultimate_coverage_table)

# =============================================================================
# TABLE 9: Functional Coverage Table (EXD)
# =============================================================================

functional_coverage_table <- data.frame(s = 1:9)
for (k in 1:9) {
    k_char <- as.character(k)
    if (k_char %in% names(all_results_sens)) {
        cov_vals <- safe_extract(all_results_sens[[k_char]]$coverage, "functional_coverage_pct", "EXD")
        functional_coverage_table[[paste0("K=", k)]] <- ifelse(
            is.na(cov_vals), "-",
            as.character(round(cov_vals, 4))
        )
    } else {
        functional_coverage_table[[paste0("K=", k)]] <- "-"
    }
}
cat("\nFunctional Coverage (EXD) by K:\n")
print(functional_coverage_table)

# =============================================================================
# TABLE 10: Ultimate Interval Score Table (EXD)
# =============================================================================

ultimate_score_table <- data.frame(s = 1:9)
for (k in 1:9) {
    k_char <- as.character(k)
    if (k_char %in% names(all_results_sens)) {
        score_vals <- safe_extract(all_results_sens[[k_char]]$coverage, "ultimate_score", "EXD")
        ultimate_score_table[[paste0("K=", k)]] <- ifelse(
            is.na(score_vals), "-",
            as.character(round(score_vals, 6))
        )
    } else {
        ultimate_score_table[[paste0("K=", k)]] <- "-"
    }
}
cat("\nUltimate Interval Score (EXD) by K:\n")
print(ultimate_score_table)

# =============================================================================
# TABLE 11: Functional Interval Score Table (EXD)
# =============================================================================

functional_score_table <- data.frame(s = 1:9)
for (k in 1:9) {
    k_char <- as.character(k)
    if (k_char %in% names(all_results_sens)) {
        score_vals <- safe_extract(all_results_sens[[k_char]]$coverage, "functional_score", "EXD")
        functional_score_table[[paste0("K=", k)]] <- ifelse(
            is.na(score_vals), "-",
            as.character(round(score_vals, 6))
        )
    } else {
        functional_score_table[[paste0("K=", k)]] <- "-"
    }
}
cat("\nFunctional Interval Score (EXD) by K:\n")
print(functional_score_table)

# =============================================================================
# FIGURE 5: CLR Tracking Plots (from fix_origin_backtest.R)
# =============================================================================
# Shows how ultimate cumulative loss predictions and intervals evolve
# as more development lags become known

cat("\n=============================================================================\n")
cat("FIGURE 5: CLR Tracking Plots for Representative Companies\n")
cat("=============================================================================\n")

# Set up variables for fix_origin plotting functions (from forecast_plot.R)
# These global variables are used by the _fix plotting functions
# NOTE: Must use the _fix versions which contain fixed-origin backtest data
test_cov_all_fix <- test_cov_3co_fix
test_matrix_all_fix <- test_matrix_3co_fix
pred_all_fix <- pred_3co_fix
boot_pred_b_all_fix <- boot_pred_b_3co_fix

# Figure 5: Three tracking plots side by side
fig5_P53 <- predict_clr_tracking_fix("P53")
fig5_C260 <- predict_clr_tracking_fix("C260")
fig5_P1406 <- predict_clr_tracking_fix("P1406")

# Arrange in a row
grid.arrange(fig5_P53, fig5_C260, fig5_P1406, ncol = 3)

# =============================================================================
# FIGURE 9: Grid Plots (ILR and CLR) for Representative Companies
# =============================================================================
# Shows prediction plots for all lags m=1 to 9
# Left column: Incremental Loss Ratios (ILR)
# Right column: Cumulative Loss Ratios (CLR)

cat("\n=============================================================================\n")
cat("FIGURE 9: Grid Plots (ILR and CLR) for Representative Companies\n")
cat("=============================================================================\n")

# Figure 9(a): Company P53 ILR
fig9a_P53_ilr <- predict_plot_grid_fix("P53", type = "i", lb = NULL, ub = NULL)
# Figure 9(b): Company P53 CLR (auto-scale since ep is normalized to 1.0)
fig9b_P53_clr <- predict_plot_grid_fix("P53", type = "c", lb = NULL, ub = NULL)

# Figure 9(c): Company C260 ILR
fig9c_C260_ilr <- predict_plot_grid_fix("C260", type = "i", lb = NULL, ub = NULL)
# Figure 9(d): Company C260 CLR (auto-scale since ep is normalized to 1.0)
fig9d_C260_clr <- predict_plot_grid_fix("C260", type = "c", lb = NULL, ub = NULL)

# Figure 9(e): Company P1406 ILR
fig9e_P1406_ilr <- predict_plot_grid_fix("P1406", type = "i", lb = NULL, ub = NULL)
# Figure 9(f): Company P1406 CLR (auto-scale since ep is normalized to 1.0)
fig9f_P1406_clr <- predict_plot_grid_fix("P1406", type = "c", lb = NULL, ub = NULL)

# Display Figure 9 plots
cat("  Figure 9(a): Company P53 ILR\n")
print(fig9a_P53_ilr)
cat("  Figure 9(b): Company P53 CLR\n")
print(fig9b_P53_clr)
cat("  Figure 9(c): Company C260 ILR\n")
print(fig9c_C260_ilr)
cat("  Figure 9(d): Company C260 CLR\n")
print(fig9d_C260_clr)
cat("  Figure 9(e): Company P1406 ILR\n")
print(fig9e_P1406_ilr)
cat("  Figure 9(f): Company P1406 CLR\n")
print(fig9f_P1406_clr)

# =============================================================================
# TABLE 12: Residual Function i.i.d. Test Results
# =============================================================================

cat("\n=============================================================================\n")
cat("TABLE 12: Residual Function i.i.d. Test Results\n")
cat("=============================================================================\n")

table12 <- data.frame(
    m = eb_test_results_fix$summary$m,
    K_star = eb_test_results_fix$summary$K,
    n = eb_test_results_fix$summary$n_obs,
    Shapiro = paste0(eb_test_results_fix$summary$shapiro_reject, "/10"),
    Ljung_Box = paste0(eb_test_results_fix$summary$ljung_reject, "/10"),
    Runs = paste0(eb_test_results_fix$summary$runs_reject, "/10"),
    Avg_Cross_Cor = round(eb_test_results_fix$summary$avg_cross_cor, 3)
)
colnames(table12) <- c("m", "K*", "n", "Shapiro", "Ljung-Box", "Runs", "Avg. Cross-Cor")
print(table12)

# =============================================================================
# TABLE 6: Optimal number of principal components and regularization parameters
# =============================================================================

cat("\n=============================================================================\n")
cat("TABLE 6: Optimal Parameters (Model for Full Analysis)\n")
cat("=============================================================================\n")

optim_param_table <- data.frame(
    ay_dev_lag = paste0(2010:2002, " (", 1:9, ")"),
    n_train = n_train_all,
    K = k_all,
    lambda = round(lam_all, 4),
    MAPE = round(mape_all, 4)
)
print(optim_param_table)

# =============================================================================
# FIGURE 6: Prediction Plots for Example Companies
# =============================================================================
# Note: These use the 3-company subset data with normalized EP
# Y-axis shows loss ratios (not dollar amounts) since EP is normalized

cat("\n=============================================================================\n")
cat("FIGURE 6: Generating prediction plots for example companies...\n")
cat("=============================================================================\n")

# Set up plot variables for predict_plot function
# These need to reference the 3-company data
test_cov_all <- test_cov_3co
test_matrix_all <- test_matrix_3co
pred_all <- pred_3co
boot_pred_b_all <- boot_pred_b_3co

cat("Note: Plots use loss ratios (normalized EP = 1.0)\n")
cat("  Figure 6(a): P53, m=3\n")
cat("  Figure 6(b): C260, m=4\n")
cat("  Figure 6(c): P1406, m=6\n")

# Generate the plots
predict_plot(3, "P53")$plt_c
predict_plot(4, "C260")$plt_c
predict_plot(6, "P1406")$plt_c

# =============================================================================
# TABLE 7: Ultimate Scoring Metrics
# =============================================================================
# Pre-computed in replication_results.Rdata as ult_metrics_summary

cat("\n=============================================================================\n")
cat("TABLE 7: Ultimate Scoring Metrics\n")
cat("=============================================================================\n")

if (exists("ult_metrics_summary")) {
    table7 <- data.frame(
        Lag_s = 1:9,
        MAPE_PLS = sprintf("%.2f%%", ult_metrics_summary[, "PLS_MAPE"] * 100),
        MAPE_CL = sprintf("%.2f%%", ult_metrics_summary[, "CL_MAPE"] * 100),
        Cov_PLS = sprintf("%.2f%%", ult_metrics_summary[, "PLS_cov"] * 100),
        Cov_EXD = sprintf("%.2f%%", ult_metrics_summary[, "EXD_cov"] * 100),
        Cov_CL = sprintf("%.2f%%", ult_metrics_summary[, "CL_cov"] * 100),
        IS_PLS = round(ult_metrics_summary[, "PLS_score"], 4),
        IS_EXD = round(ult_metrics_summary[, "EXD_score"], 4),
        IS_CL = round(ult_metrics_summary[, "CL_score"], 4)
    )
    colnames(table7) <- c(
        "Lag s", "MAPE PLS", "MAPE CL", "Cov% PLS", "Cov% EXD", "Cov% CL",
        "IS PLS", "IS EXD", "IS CL"
    )
    print(table7)
} else {
    cat("Note: ult_metrics_summary not found in replication_results.Rdata\n")
    cat("This table requires pre-computed ultimate metrics.\n")
}

# =============================================================================
# TABLE 8: Functional Scoring Metrics
# =============================================================================

cat("\n=============================================================================\n")
cat("TABLE 8: Functional Scoring Metrics\n")
cat("=============================================================================\n")

if (exists("ult_metrics_summary")) {
    table8 <- data.frame(
        Lag_s = 1:9,
        Cov_PLS = sprintf("%.2f%%", ult_metrics_summary[, "PLS_func_cov"] * 100),
        Cov_EXD = sprintf("%.2f%%", ult_metrics_summary[, "EXD_func_cov"] * 100),
        Cov_CL = sprintf("%.2f%%", ult_metrics_summary[, "CL_func_cov"] * 100),
        IS_PLS = round(ult_metrics_summary[, "PLS_func_score"], 4),
        IS_EXD = round(ult_metrics_summary[, "EXD_func_score"], 4),
        IS_CL = round(ult_metrics_summary[, "CL_func_score"], 4)
    )
    colnames(table8) <- c("s", "Cov% PLS", "Cov% EXD", "Cov% CL", "IS PLS", "IS EXD", "IS CL")
    print(table8)
} else {
    cat("Note: ult_metrics_summary not found in replication_results.Rdata\n")
    cat("This table requires pre-computed functional metrics.\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=============================================================================\n")
cat("REPLICATION COMPLETE\n")
cat("=============================================================================\n")
cat("Tables generated:\n")
cat("  - Table 7: MAPE by K value\n")
cat("  - Table 8: Ultimate Coverage (EXD) by K\n")
cat("  - Table 9: Functional Coverage (EXD) by K\n")
cat("  - Table 10: Ultimate Interval Score (EXD) by K\n")
cat("  - Table 11: Functional Interval Score (EXD) by K\n")
cat("  - Table 12: Residual i.i.d. tests\n")
cat("  - Table 6: Optimal parameters\n")
cat("  - Table 7: Ultimate Scoring Metrics (PLS vs CL)\n")
cat("  - Table 8: Functional Scoring Metrics (PLS vs CL)\n")
cat("\nFigures generated:\n")
cat("  - Figure 5: CLR Tracking plots (P53, C260, P1406)\n")
cat("  - Figure 6: Prediction plots (P53, C260, P1406)\n")
cat("  - Figure 9: Grid plots with ILR and CLR (P53, C260, P1406)\n")
cat("=============================================================================\n")
