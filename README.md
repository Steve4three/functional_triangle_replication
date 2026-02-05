# Functional Triangle Replication Package

This repository contains the replication code for the Functional Principal Component Analysis (FPCA) approach to loss reserving in property-casualty insurance.

## Requirements

### R Version
- R 4.0.0 or higher

### Required Packages
```r
install.packages(c("tidyverse", "gridExtra", "fdaoutlier", "glmnet", "caret", "psych", "parallel", "doParallel"))
```

## Quick Start: Running the Replication

To reproduce all tables and figures from the paper:

```r
source("replication.R")
```

This script uses pre-computed results and does not require access to the raw proprietary data.

## Expected Output

Running `replication.R` will generate:

### Tables
| Table | Description |
|-------|-------------|
| Table 6 | Optimal parameters (K, λ, MAPE) for full analysis by development lag |
| Table 7 | MAPE comparison across all K values (1-9) for sensitivity analysis |
| Table 8 | Ultimate coverage by K value (EXD method) |
| Table 9 | Functional coverage by K value (EXD method) |
| Table 10 | Ultimate interval score by K value (EXD method) |
| Table 11 | Functional interval score by K value (EXD method) |
| Table 12 | Residual function i.i.d. test results (Shapiro, Ljung-Box, Runs) |

### Figures
| Figure | Description |
|--------|-------------|
| Figure 5 | CLR tracking plots showing ultimate cumulative loss predictions evolving as more lags become known |
| Figure 6 | Prediction plots with 95% confidence intervals for representative companies |
| Figure 9 | Grid plots showing ILR (Incremental Loss Ratio) and CLR (Cumulative Loss Ratio) for each lag m=1-9 |

**Note:** Plots use loss ratios (not dollar amounts) since EP is normalized to 1.0 in the replication data.

## File Structure

### Replication Files (Use These)
| File | Description |
|------|-------------|
| `replication.R` | Main replication script - run this to reproduce all results |
| `replication_data.Rdata` | Pre-computed data for 3 representative companies (P53, C260, P1406) |
| `replication_results.Rdata` | Pre-computed aggregate results for tables |
| `sensitivity_K_results.Rdata` | Sensitivity analysis results for K values 1-9 |

### Function Files (Methodology Implementation)

These files contain the methodology implementation and can be reviewed for understanding the approach:

#### Core Functions
| File | Description | Requires Raw Data |
|------|-------------|-------------------|
| `functions.R` | Core FDA functions: `data_gen()`, `data_split()`, `fold_gen()`, `beta_reg()`, `pls_triangle()` | No (sourced by replication.R) |
| `tuning.R` | Hyperparameter tuning: `k_tune()` for K selection, `l_tune()` for λ selection | Yes (runs on full dataset) |
| `chain_ladder.R` | Chain ladder baseline: `cl()` for traditional method, `calculate_reserves()` | Yes (runs on full dataset) |
| `prediction_interval.R` | Interval functions: `pls_interval()`, `depth_interval()`, `depth_interval_fix()` | No (sourced by replication.R) |
| `evaluation.R` | Evaluation metrics: `calculate_coverage_score()`, `calculate_ultimate_metrics()` | No (sourced by replication.R) |
| `forecast_plot.R` | Visualization: `predict_plot()`, `predict_plot_grid_fix()`, `predict_clr_tracking_fix()` | No (sourced by replication.R) |

#### Data Processing (Require Raw Data)
| File | Description |
|------|-------------|
| `data_prep.R` | Data loading and preprocessing: `read_data()`, `prepare_data()`, `create_pca_matrix()` |
| `eda_outliers.R` | Outlier detection: `devlag_summary()`, `plot_bag_outliers()`, functional boxplots |
| `eda_covariates.R` | EDA visualizations: `plot_bag_covariate()`, `plot_median_accident_years()` |
| `model.R` | Full model pipeline: `process_m()` for running the complete analysis |
| `fix_origin_backtest.R` | Fixed-origin backtest: `process_m_fix()`, `k_tune_fix()`, `l_tune_fix()` |
| `sensitivity_K.R` | K sensitivity analysis: `process_m_sens()`, `l_tune_sens()` |

### Function Details

#### `functions.R` - Core FDA Implementation

```r
# Generate training/test splits by accident year cutoffs
data_gen(data, cutoff1, cutoff2)

# Random train/test split for cross-validation
data_split(data, p)

# K-fold cross-validation data generation
fold_gen(data, fn = 5)

# Beta regression using LASSO or stepwise selection
beta_reg(beta_data, covariate, k, alpha = NULL, lasso = TRUE)

# PLS triangle prediction combining FDA with PCA decomposition
pls_triangle(lambda, m, beta_model, covariate, test_matrix, pca, k, lasso = TRUE)
```

#### `tuning.R` - Hyperparameter Selection

```r
# Select optimal K (number of principal components) via 5-fold CV
k_tune(m, data, n_fold = 5)
# Returns: list(mean = vector of MAPE by K, msd = mean + sd for one-SE rule)

# Tune regularization parameter λ given fixed K
l_tune(m, k, data, n_fold = 5)
# Returns: list(lam = optimal λ, mape = resulting MAPE)
```

#### `prediction_interval.R` - Uncertainty Quantification

```r
# Pointwise prediction intervals using bootstrap percentiles
pls_interval(m, company, test_cov_list, test_matrix_list, pred_list, boot_pred_b_list)

# Functional prediction intervals using depth-based methods
depth_interval(m, company, test_cov_list, test_matrix_list, pred_list, boot_pred_b_list, method = "extremal", ef = 2.5)

# Fixed-origin backtest version (uses global _fix variables)
depth_interval_fix(m, company, method = "extremal", ef = 2.5)
```

#### `evaluation.R` - Performance Assessment

```r
# Calculate coverage and interval scores for all companies
calculate_coverage_score(company_list, test_cov_list, test_matrix_list, pred_list, boot_pred_b_list, valid_m = rep(TRUE, 9))

# Detailed metrics including comparison with Chain Ladder
calculate_ultimate_metrics(company_list, test_cov_list, test_matrix_list, pred_list, boot_pred_b_list)
```

### How the Non-Replicable Scripts Work

These scripts require access to the full raw dataset which is not shared due to confidentiality:

#### `data_prep.R`
Loads and processes raw SNL data:
1. Reads triangle data (paid/reported/incurred losses by development lag)
2. Joins with company characteristics
3. Calculates incremental loss ratios from cumulative data
4. Standardizes categorical variables (business focus, ownership, geography)

#### `model.R`
Runs the complete analysis pipeline for all companies:
1. Loops through each development lag m (1-9)
2. Tunes K and λ for each lag
3. Generates point predictions and bootstrap samples
4. Updates imputed data progressively

#### `fix_origin_backtest.R`
Performs fixed-origin validation (training on pre-2010 data, testing on accident year 2010):
1. Trains models with progressively revealed lags
2. Evaluates prediction accuracy as more information becomes available
3. Generates results used in the sensitivity tables

#### `sensitivity_K.R`
Analyzes model stability across different K choices:
1. Fixes K at values 1-9 (instead of tuning)
2. Tunes only λ for each fixed K
3. Computes coverage and scores to assess robustness

## Data Privacy

The shared data files contain only:
- **3 anonymized company identifiers** (P53, C260, P1406)
- **Normalized earned premiums** (EP = 1.0 for all observations)
- **Loss ratios** (between 0 and 1, not dollar amounts)
- **Aggregate model statistics** (K, λ, MAPE, coverage rates)

No raw loss dollar amounts or identifiable company information is included.

## Citation

If you use this code, please cite: [Paper citation to be added]

## Contact

For questions about the methodology or code, please contact: [Contact information to be added]

## License

[License to be added]
