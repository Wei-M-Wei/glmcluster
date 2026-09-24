rm(list = ls())

# =============================================================================
# PACKAGES AND FUNCTIONS
# =============================================================================

library(sandwich)
library(gtools)
library(foreach)
library(parallel)
library(doParallel)
library(doRNG)
library(pcluster)
library(fixest)
library(alpaca)
library(dplyr)
library(NNR)

source("functions for nonlinear.R")


# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

sample_file <- "Sample_use.csv"

# Column layout:
#   1: id
#   2: time
#   3: outcome
#   4+: covariates

header_names <- names(
  read.csv(sample_file, nrows = 0, check.names = FALSE)
)

col_classes <- rep(NA, length(header_names))
col_classes[1:2] <- "character"

df <- read.csv(
  sample_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = col_classes
)

id_var      <- names(df)[1]
time_var    <- names(df)[2]
outcome_var <- names(df)[3]
covs        <- names(df)[4:ncol(df)]

cat("Dataset dimensions:", dim(df), "\n")
cat("ID variable:", id_var, "\n")
cat("Time variable:", time_var, "\n")
cat("Outcome variable:", outcome_var, "\n")
cat("Number of ids:", length(unique(df[[id_var]])), "\n")
cat("Number of time periods:", length(unique(df[[time_var]])), "\n")
cat("Time range:", min(df[[time_var]]), "-", max(df[[time_var]]), "\n")
cat("Variable names:\n")
print(names(df))
cat("\n")


# =============================================================================
# 2. DATA CHECKS
# =============================================================================

cat("Unique time periods:\n")
print(sort(unique(df[[time_var]])))

cat("\nUnique ids (first 10):\n")
print(head(sort(unique(df[[id_var]])), 10))

# Missing values
missing_cols <- names(df)[colSums(is.na(df)) > 0]

if (length(missing_cols) > 0) {
  stop(
    "Missing values found in: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Duplicate id-time observations
duplicate_rows <- duplicated(df[, c(id_var, time_var)])

if (any(duplicate_rows)) {
  stop(
    "Duplicate id-time rows found: ",
    sum(duplicate_rows)
  )
}

# No filtering
df_final <- df

cat("\nFinal dataset dimensions:", dim(df_final), "\n")
cat(
  "Final number of ids:",
  length(unique(df_final[[id_var]])),
  "\n\n"
)


# =============================================================================
# 3. CREATE UNIT AND TIME INDICES
# =============================================================================

df_vars <- df_final[
  ,
  c(id_var, time_var, outcome_var, covs),
  drop = FALSE
]

time_values <- unique(df_vars[[time_var]])
time_values_numeric <- suppressWarnings(as.numeric(time_values))

time_levels <- if (all(!is.na(time_values_numeric))) {
  time_values[order(time_values_numeric)]
} else {
  sort(time_values)
}

df_vars$unit_id <- as.integer(
  factor(
    df_vars[[id_var]],
    levels = unique(df_vars[[id_var]])
  )
)

df_vars$time_id <- as.integer(
  factor(
    df_vars[[time_var]],
    levels = time_levels
  )
)

cat(
  "Covariates included:",
  paste(covs, collapse = ", "),
  "\n\n"
)


# =============================================================================
# 4. VERIFY BALANCED PANEL
# =============================================================================

panel_N <- length(unique(df_vars$unit_id))
panel_T <- length(unique(df_vars$time_id))
expected_rows <- panel_N * panel_T

cat(
  "Panel dimensions: N =", panel_N,
  ", T =", panel_T, "\n"
)
cat("Observed rows:", nrow(df_vars), "\n")
cat("Expected rows:", expected_rows, "\n")

if (nrow(df_vars) != expected_rows) {
  stop(
    "Sample_use.csv is not balanced. Expected ",
    expected_rows,
    " rows but found ",
    nrow(df_vars),
    "."
  )
}

counts_by_id <- table(df_vars$unit_id)

if (!all(counts_by_id == panel_T)) {
  stop(
    "Not every id has exactly ",
    panel_T,
    " time periods."
  )
}

for (i in unique(df_vars$unit_id)[1:min(3, panel_N)]) {
  raw_id <- df_vars[[id_var]][
    which(df_vars$unit_id == i)[1]
  ]

  cat(
    "ID", raw_id, "time periods:",
    paste(
      df_vars[df_vars$unit_id == i, time_var],
      collapse = ", "
    ),
    "\n"
  )
}


# =============================================================================
# 5. NORMALIZE COVARIATES
# =============================================================================

df_ready <- df_vars[
  ,
  c("unit_id", "time_id", outcome_var, covs),
  drop = FALSE
]

# Keep an unnormalized copy
df_ready_orig <- df_ready

num_covs_all <- setdiff(
  names(df_ready),
  c("unit_id", "time_id", outcome_var)
)

num_covs_all <- num_covs_all[
  sapply(df_ready[num_covs_all], is.numeric)
]

# Identify dummy variables
is_dummy <- sapply(
  num_covs_all,
  function(var) {
    unique_vals <- unique(na.omit(df_ready[[var]]))

    length(unique_vals) <= 2 &&
      all(unique_vals %in% c(0, 1))
  }
)

dummy_vars <- names(is_dummy[is_dummy])
continuous_vars <- setdiff(num_covs_all, dummy_vars)

cat(
  "Dummy variables (not normalized):",
  paste(dummy_vars, collapse = ", "),
  "\n"
)

cat(
  "Continuous variables:",
  paste(continuous_vars, collapse = ", "),
  "\n\n"
)

# Standardize continuous covariates
if (length(continuous_vars) > 0) {

  df_ready[continuous_vars] <- as.data.frame(
    scale(df_ready[continuous_vars])
  )

  message(
    "Normalized continuous covariates: ",
    paste(continuous_vars, collapse = ", ")
  )

} else {

  message("No continuous covariates found to normalize.")
}

num_covs <- c(continuous_vars, dummy_vars)

df_ready <- df_ready[
  ,
  c("unit_id", "time_id", outcome_var, num_covs),
  drop = FALSE
]

cat(
  "\nAll covariates for estimation:",
  paste(num_covs, collapse = ", "),
  "\n"
)

cat(
  "Covariate count:",
  length(num_covs),
  "\n\n"
)


# =============================================================================
# 6. CONVERT DATA TO N x T MATRICES
# =============================================================================

cat("=== CONVERTING TO MATRIX FORM (N x T) ===\n")

panel_N <- length(unique(df_ready$unit_id))
panel_T <- length(unique(df_ready$time_id))

cat(
  "Panel dimensions: N =", panel_N,
  ", T =", panel_T, "\n"
)

cat(
  "Total observations:",
  panel_N * panel_T,
  "\n\n"
)

# Outcome matrix
Y_mat <- matrix(
  NA_real_,
  nrow = panel_N,
  ncol = panel_T
)

for (r in seq_len(nrow(df_ready))) {
  i <- df_ready$unit_id[r]
  t <- df_ready$time_id[r]

  Y_mat[i, t] <- df_ready[[outcome_var]][r]
}

cat(
  "Outcome matrix Y_mat:",
  nrow(Y_mat), "x", ncol(Y_mat), "\n"
)

cat(
  "No missing values:",
  all(!is.na(Y_mat)),
  "\n\n"
)

# Covariate matrices
X_list <- lapply(
  num_covs,
  function(varname) {

    X <- matrix(
      NA_real_,
      nrow = panel_N,
      ncol = panel_T
    )

    for (r in seq_len(nrow(df_ready))) {
      i <- df_ready$unit_id[r]
      t <- df_ready$time_id[r]

      X[i, t] <- df_ready[[varname]][r]
    }

    X
  }
)

names(X_list) <- num_covs

cat(
  "Created", length(X_list),
  "covariate matrices:\n"
)

for (k in seq_along(X_list)) {

  var_name <- names(X_list)[k]

  cat(
    "  -", var_name,
    "(", nrow(X_list[[k]]),
    "x", ncol(X_list[[k]]), ")",
    "- Complete:",
    all(!is.na(X_list[[k]])),
    "\n"
  )
}

cat("\nMatrix conversion complete.\n\n")


# =============================================================================
# 7. ESTIMATION SETTINGS
# =============================================================================

set.seed(1)

Y <- Y_mat
N <- nrow(Y)
T <- ncol(Y)

link <- "logit"
discrete_type <- "additive"
gamma <- 1

theta <- rep(0, length(num_covs))


# =============================================================================
# 8. BC-E: ELBOW METHOD
# =============================================================================

model_elbow <- GFE_est(
  Y_mat,
  X_list = X_list,
  link = link,
  discrete_type = discrete_type,
  gamma = gamma,
  cluster_type = "kcenter",
  group_selection = "elbow",
  unit_max_groups = floor(3 * N^(1 / 3)),
  time_max_groups = floor(3 * T^(1 / 3))
)

G_elbow <- model_elbow$G
C_elbow <- model_elbow$C

beta_elbow <- coef(model_elbow$model)
se_elbow <- model_elbow$se

beta_BC_E <- coef(model_elbow$est_analytical)
se_BC_E <- model_elbow$se_ana

cat("\n=== BC-E ===\n")
cat("G =", G_elbow, ", C =", C_elbow, "\n")

print(beta_BC_E)
print(se_BC_E)
print(beta_BC_E / se_BC_E)

cover_BC_E <-
  (beta_BC_E - 1.645 * se_BC_E <= theta) *
  (beta_BC_E + 1.645 * se_BC_E >= theta)

print(cover_BC_E)


# =============================================================================
# 9. BC-N: NOISE-FLOOR RULE
# =============================================================================

model_noise <- GFE_est(
  Y_mat,
  X_list = X_list,
  link = link,
  discrete_type = discrete_type,
  gamma = gamma,
  cluster_type = "kcenter",
  unit_max_groups = floor(3 * N^(1 / 3)),
  time_max_groups = floor(3 * T^(1 / 3))
)

G_noise <- model_noise$G
C_noise <- model_noise$C

beta_noise <- coef(model_noise$model)
se_noise <- model_noise$se

beta_BC_N <- coef(model_noise$est_analytical)
se_BC_N <- model_noise$se_ana

cat("\n=== BC-N ===\n")
cat("G =", G_noise, ", C =", C_noise, "\n")

print(beta_BC_N)
print(se_BC_N)
print(beta_BC_N / se_BC_N)

cover_BC_N <-
  (beta_BC_N - 1.645 * se_BC_N <= theta) *
  (beta_BC_N + 1.645 * se_BC_N >= theta)

print(cover_BC_N)


# =============================================================================
# 10. CROSS-FITTING: ELBOW METHOD
# =============================================================================

model_cf_elbow <- GFE_est_cf(
  Y_mat,
  X_list = X_list,
  link = link,
  discrete_type = discrete_type,
  gamma = gamma,
  cluster_type = "kcenter",
  group_selection = "elbow",
  unit_max_groups = 3 * floor(N^(1 / 3)),
  time_max_groups = 3 * floor(T^(1 / 3))
)

# Pooled cross-fitting estimator
beta_BC_CFP_E <- coef(model_cf_elbow$est_analytical)
se_BC_CFP_E <- model_cf_elbow$se_ana

# Fold-averaged cross-fitting estimator
beta_BC_CFA_E <- model_cf_elbow$beta_fold_avg
se_BC_CFA_E <- model_cf_elbow$se_fold_avg

cat("\n=== Cross-fitting: elbow ===\n")
cat(
  "G =", model_cf_elbow$G,
  ", C =", model_cf_elbow$C,
  "\n"
)

cat("\nBC-CFP-E:\n")
print(beta_BC_CFP_E)
print(se_BC_CFP_E)

cat("\nBC-CFA-E:\n")
print(beta_BC_CFA_E)
print(se_BC_CFA_E)


# =============================================================================
# 11. CROSS-FITTING: NOISE-FLOOR RULE
# =============================================================================

model_cf_noise <- GFE_est_cf(
  Y_mat,
  X_list = X_list,
  link = link,
  discrete_type = discrete_type,
  gamma = gamma,
  cluster_type = "kcenter",
  unit_max_groups = 3 * floor(N^(1 / 3)),
  time_max_groups = 3 * floor(T^(1 / 3))
)

# Pooled cross-fitting estimator
beta_BC_CFP_N <- coef(model_cf_noise$est_analytical)
se_BC_CFP_N <- model_cf_noise$se_ana

# Fold-averaged cross-fitting estimator
beta_BC_CFA_N <- model_cf_noise$beta_fold_avg
se_BC_CFA_N <- model_cf_noise$se_fold_avg

cat("\n=== Cross-fitting: noise-floor ===\n")
cat(
  "G =", model_cf_noise$G,
  ", C =", model_cf_noise$C,
  "\n"
)

cat("\nBC-CFP-N:\n")
print(beta_BC_CFP_N)
print(se_BC_CFP_N)

cat("\nBC-CFA-N:\n")
print(beta_BC_CFA_N)
print(se_BC_CFA_N)


# =============================================================================
# 12. IFE ANALYTICAL BIAS CORRECTION
# =============================================================================

est_IFE <- NNRPanel_estimate(
  data_frame = as.matrix(df_ready),
  func = "logit",
  R_max = 3,
  delta = 1.05,
  iter_max = 10000,
  tol = 1e-8,
  delta_1 = 0.5
)

beta_IFE_ABC <- est_IFE$beta_corr_data
se_IFE_ABC <- est_IFE$std_beta_raw_data

cat("\n=== IFE-ABC ===\n")
print(beta_IFE_ABC)
print(se_IFE_ABC)
print(beta_IFE_ABC / se_IFE_ABC)


# =============================================================================
# 13. TWFE ANALYTICAL BIAS CORRECTION
# =============================================================================

x_vars <- colnames(df_ready)[4:ncol(df_ready)]

fml_twfe <- as.formula(
  paste(
    outcome_var,
    "~",
    paste(x_vars, collapse = " + "),
    "| unit_id + time_id"
  )
)

est_TWFE <- alpaca::feglm(
  fml_twfe,
  data = data.frame(df_ready),
  family = binomial(link = "logit"),
  control = feglmControl(
    dev.tol = 1e-6,
    center.tol = 1e-6,
    iter.max = 1000,
    drop.pc = FALSE
  )
)

BC_TWFE <- biasCorr(
  est_TWFE,
  L = 0,
  panel.structure = "classic"
)

beta_TWFE_ABC <- coef(BC_TWFE)
se_TWFE_ABC <- sqrt(diag(solve(BC_TWFE$Hessian)))

cat("\n=== TWFE-ABC ===\n")
print(beta_TWFE_ABC)
print(se_TWFE_ABC)
print(beta_TWFE_ABC / se_TWFE_ABC)


# =============================================================================
# 14. SUMMARY
# =============================================================================

results <- list(
  BC_E = list(
    estimate = beta_BC_E,
    se = se_BC_E
  ),
  BC_CFA_E = list(
    estimate = beta_BC_CFA_E,
    se = se_BC_CFA_E
  ),
  BC_CFP_E = list(
    estimate = beta_BC_CFP_E,
    se = se_BC_CFP_E
  ),
  BC_N = list(
    estimate = beta_BC_N,
    se = se_BC_N
  ),
  BC_CFA_N = list(
    estimate = beta_BC_CFA_N,
    se = se_BC_CFA_N
  ),
  BC_CFP_N = list(
    estimate = beta_BC_CFP_N,
    se = se_BC_CFP_N
  ),
  IFE_ABC = list(
    estimate = beta_IFE_ABC,
    se = se_IFE_ABC
  ),
  TWFE_ABC = list(
    estimate = beta_TWFE_ABC,
    se = se_TWFE_ABC
  )
)

results
