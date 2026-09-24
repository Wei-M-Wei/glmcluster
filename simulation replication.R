rm(list = ls())

library(sandwich)
library(gtools)
library(foreach)
library(parallel)
library(doParallel)
library(fixest)
library(alpaca)
library(NNR)

# Functions

source("DGP truncated.R")
source("functions for nonlinear.R")

# Parallel setup

# On RStudio Server, use the number of processes allocated by SLURM.
# On a local PC, fall back to detectCores() - 1.
slurm_workers <- suppressWarnings(as.integer(Sys.getenv("SLURM_NTASKS")))

if (!is.na(slurm_workers) && slurm_workers > 0) {
  nworkers <- slurm_workers
} else {
  nworkers <- max(1, detectCores() - 1)
}

cat("Number of parallel workers:", nworkers, "\n")

cl <- makeCluster(50)
registerDoParallel(cl)

# Force each Monte Carlo worker to use one BLAS/OpenMP thread
clusterEvalQ(cl, {

  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  }

  NULL
})

# Packages required by workers
packages_to_export <- c("fixest", "alpaca", "sandwich", "NNR")

# Parameters

N <- 50
T <- 20

dimtheta <- 1
K <- dimtheta

theta <- -0.5

rho <- 0.5
kappa <- 0.5

Rep <- 1000

DGP <- 1
nu <- 5

link <- "probit"

cluster_type <- "kcenter"

# IFE estimates with absolute value above this bound are treated as numerical failures
IFE_bound <- 100

# CF fold-average estimates with absolute value above this bound are treated as numerical failures
CF_foldavg_bound <- 100

# Helper functions

first_value <- function(x) {

  x <- c(x)

  if (length(x) == 0) {
    NA_real_
  } else {
    x[1]
  }
}

cover_first <- function(beta, se, theta0) {

  if (!is.finite(beta) || !is.finite(se)) {
    return(NA_real_)
  }

  as.numeric(beta - 1.96 * se <= theta0 && beta + 1.96 * se >= theta0)
}

# Simulation

for (theta0 in c(-0.5)) {

  for (link in c("logit", "probit")) {

    for (DGP in c(1,2,3,4)) {

      for (T in c(20, 30, 40, 50)) {

        # Parameter dimension

        if (DGP == 4) {

          dimtheta <- 2
          K <- dimtheta

          theta <- c(theta0, theta0)

        } else {

          dimtheta <- 1
          K <- dimtheta

          theta <- theta0
        }

        # We evaluate the first coefficient
        target_theta <- theta[1]

        cat("\nStarting:", "N =", N, "T =", T, "DGP =", DGP, "link =", link, "theta =", paste(theta, collapse = ", "), "\n")

        # Monte Carlo

        start_time <- Sys.time()

        result <- foreach(
          k = 1:Rep,
          .errorhandling = "remove",
          .packages = packages_to_export
        ) %dopar% {

          set.seed(k)

          # Generate data

          data <- generate(N = N, T = T, DGP = DGP, rho = rho, kappa = kappa, theta = theta, link = link, nu = nu)

          Y <- reshape_to_matrix(data, "id", "time", "vY")

          x_vars <- if (
            all(c("vX1", "vX2") %in% names(data))
          ) {
            c("vX1", "vX2")
          } else {
            "vX"
          }

          X_list <- lapply(
            x_vars,
            function(var) {
              reshape_to_matrix(data, "id", "time", var)
            }
          )

          # 1. Interaction GFE

          model <- GFE_est(Y, X_list, link = link, discrete_type = "interaction", gamma = 1, cluster_type = cluster_type, compute_analytical = FALSE)

          G <- model$G
          C <- model$C

          cof <- first_value(coef(model$model))

          se_interaction <- first_value(model$se)

          cover_interaction <- cover_first(cof, se_interaction, target_theta)

          # 2. Additive GFE

          model_additive <- GFE_est(Y, X_list = X_list, link = link, discrete_type = "additive", gamma = 1, cluster_type = cluster_type, reuse_analytical_fit = TRUE, unit_max_groups = 3 * floor(N^(1 / 3)), time_max_groups = 3 * floor(T^(1 / 3)))

          G_additive <- model_additive$G
          C_additive <- model_additive$C

          beta_additive <- first_value(coef(model_additive$model))

          se_additive <- first_value(model_additive$se)

          cover_additive <- cover_first(beta_additive, se_additive, target_theta)

          # 3. Additive analytical correction

          cof_corrected <- first_value(coef(model_additive$est_analytical))

          se_analytical <- first_value(model_additive$se_ana)

          cover_analytical <- cover_first(cof_corrected, se_analytical, target_theta)

          # 4. IFE analytical correction

          est_IFE <- tryCatch(
            NNRPanel_estimate(data_frame = data, func = link, delta = 1.05, R_max = 5, s = 100, iter_max = 10000, tol = 1e-6, delta_1 = 0.5),
            error = function(e) NULL
          )

          if (is.null(est_IFE)) {

            cof_IFE <- NA_real_
            se_IFE <- NA_real_
            cover_IFE <- NA_real_
            IFE_failed <- TRUE

          } else {

            cof_IFE_raw <- first_value(est_IFE$beta_corr_data)
            se_IFE_raw <- first_value(est_IFE$std_beta_corr_data)

            IFE_failed <- !is.finite(cof_IFE_raw) ||
              !is.finite(se_IFE_raw) ||
              se_IFE_raw <= 0 ||
              abs(cof_IFE_raw) > IFE_bound

            if (IFE_failed) {
              cof_IFE <- NA_real_
              se_IFE <- NA_real_
              cover_IFE <- NA_real_
            } else {
              cof_IFE <- cof_IFE_raw
              se_IFE <- se_IFE_raw
              cover_IFE <- cover_first(cof_IFE, se_IFE, target_theta)
            }
          }

          # 5. TWFE analytical correction

          fml_twfe <- as.formula(paste("vY ~", paste(x_vars, collapse = " + "), "| id + time"))

          est_TWFE <- alpaca::feglm(fml_twfe, data = data, family = binomial(link = link), control = feglmControl(dev.tol = 1e-6, center.tol = 1e-6, iter.max = 1000, drop.pc = FALSE))

          BC_TWFE <- biasCorr(est_TWFE, L = 0, panel.structure = "classic")

          cof_TWFE <- first_value(BC_TWFE$coefficients)

          se_TWFE <- first_value(sqrt(diag(solve(BC_TWFE$Hessian))))

          cover_TWFE <- cover_first(cof_TWFE, se_TWFE, target_theta)

          # 6. Cross-fitted additive estimator

          model_additive_cf <- GFE_est_cf(Y, X_list = X_list, link = link, discrete_type = "additive", gamma = 1, cluster_type = cluster_type, unit_max_groups = 3 * floor(N^(1 / 3)), time_max_groups = 3 * floor(T^(1 / 3)), compute_fold_average = TRUE, compute_pooled_uncorrected = FALSE)

          G_analytical_cf <- model_additive_cf$G
          C_analytical_cf <- model_additive_cf$C

          beta_analytical_cf <- first_value(coef(model_additive_cf$est_analytical))

          se_analytical_cf <- first_value(model_additive_cf$se_ana)

          cover_analytical_cf <- cover_first(beta_analytical_cf, se_analytical_cf, target_theta)

          # Fold-average quantities are retained
          # although they are not included in the main table.

          beta_analytical_cf_foldavg_raw <- first_value(model_additive_cf$beta_fold_avg)

          se_analytical_cf_foldavg_raw <- first_value(model_additive_cf$se_fold_avg)

          CF_foldavg_failed <- !is.finite(beta_analytical_cf_foldavg_raw) ||
            !is.finite(se_analytical_cf_foldavg_raw) ||
            se_analytical_cf_foldavg_raw <= 0 ||
            abs(beta_analytical_cf_foldavg_raw) > CF_foldavg_bound

          if (CF_foldavg_failed) {
            beta_analytical_cf_foldavg <- NA_real_
            se_analytical_cf_foldavg <- NA_real_
            cover_analytical_cf_foldavg <- NA_real_
          } else {
            beta_analytical_cf_foldavg <- beta_analytical_cf_foldavg_raw
            se_analytical_cf_foldavg <- se_analytical_cf_foldavg_raw
            cover_analytical_cf_foldavg <- cover_first(beta_analytical_cf_foldavg, se_analytical_cf_foldavg, target_theta)
          }

          # Return replication results

          results <- list(beta_interaction = cof, beta_additive = beta_additive, beta_analytical = cof_corrected, beta_analytical_cf = beta_analytical_cf, beta_analytical_cf_foldavg = beta_analytical_cf_foldavg, beta_IFE = cof_IFE, IFE_failed = IFE_failed, CF_foldavg_failed = CF_foldavg_failed, G = G, C = C, G_additive = G_additive, C_additive = C_additive, G_analytical_cf = G_analytical_cf, C_analytical_cf = C_analytical_cf, cover_interaction = cover_interaction, cover_additive = cover_additive, cover_analytical = cover_analytical, cover_analytical_cf = cover_analytical_cf, cover_analytical_cf_foldavg = cover_analytical_cf_foldavg, cover_IFE = cover_IFE, se_interaction = se_interaction, se_additive = se_additive, se_analytical = se_analytical, se_analytical_cf = se_analytical_cf, se_analytical_cf_foldavg = se_analytical_cf_foldavg, se_IFE = se_IFE, beta_TWFE = cof_TWFE, se_TWFE = se_TWFE, cover_TWFE = cover_TWFE)

          return(results)
        }

        end_time <- Sys.time()

        cat("Monte Carlo elapsed time:", round(as.numeric(difftime(end_time, start_time, units = "secs")), 2), "seconds\n")

        # Clean results

        result <- Filter(Negate(is.null), result)

        cat("Successful replications:", length(result), "out of", Rep, "\n")

        # Seven estimators:
        #
        # 1 interaction
        # 2 additive
        # 3 additive analytical correction
        # 4 analytical CF
        # 5 analytical CF fold average
        # 6 IFE
        # 7 TWFE

        estimate <- matrix(NA_real_, nrow = length(result), ncol = 7)

        cover <- matrix(NA_real_, nrow = length(result), ncol = 7)

        se <- matrix(NA_real_, nrow = length(result), ncol = 7)

        cluster_G <- matrix(NA_real_, nrow = length(result), ncol = 7)

        cluster_C <- matrix(NA_real_, nrow = length(result), ncol = 7)

        IFE_failed <- rep(NA, length(result))

        CF_foldavg_failed <- rep(NA, length(result))

        # Collect replication results

        for (i in seq_along(result)) {

          estimate[i, ] <- c(result[[i]]$beta_interaction, result[[i]]$beta_additive, result[[i]]$beta_analytical, result[[i]]$beta_analytical_cf, result[[i]]$beta_analytical_cf_foldavg, result[[i]]$beta_IFE, result[[i]]$beta_TWFE)

          cover[i, ] <- c(result[[i]]$cover_interaction, result[[i]]$cover_additive, result[[i]]$cover_analytical, result[[i]]$cover_analytical_cf, result[[i]]$cover_analytical_cf_foldavg, result[[i]]$cover_IFE, result[[i]]$cover_TWFE)

          se[i, ] <- c(result[[i]]$se_interaction, result[[i]]$se_additive, result[[i]]$se_analytical, result[[i]]$se_analytical_cf, result[[i]]$se_analytical_cf_foldavg, result[[i]]$se_IFE, result[[i]]$se_TWFE)

          cluster_G[i, ] <- c(result[[i]]$G, result[[i]]$G_additive, result[[i]]$G_additive, result[[i]]$G_analytical_cf, result[[i]]$G_analytical_cf, NA_real_, NA_real_)

          cluster_C[i, ] <- c(result[[i]]$C, result[[i]]$C_additive, result[[i]]$C_additive, result[[i]]$C_analytical_cf, result[[i]]$C_analytical_cf, NA_real_, NA_real_)

          IFE_failed[i] <- result[[i]]$IFE_failed

          CF_foldavg_failed[i] <- result[[i]]$CF_foldavg_failed
        }

        # Monte Carlo statistics

        sd_est <- apply(estimate, 2, sd, na.rm = TRUE)

        # Median estimated standard error
        se_average <- apply(se, 2, mean, na.rm = TRUE)

        # IMPORTANT:
        # all methods here are evaluated for theta[1]
        bias <- apply(estimate, 2, mean, na.rm = TRUE) - target_theta

        bias_median <- apply(estimate, 2, median, na.rm = TRUE) - target_theta

        cover_rate <- apply(cover, 2, mean, na.rm = TRUE)

        sd_ratio <- se_average / sd_est

        cluster_G_average <- c(colMeans(cluster_G[, 1:5, drop = FALSE], na.rm = TRUE), NA_real_, NA_real_)

        cluster_C_average <- c(colMeans(cluster_C[, 1:5, drop = FALSE], na.rm = TRUE), NA_real_, NA_real_)

        CF_foldavg_failure_rate <- c(NA_real_, NA_real_, NA_real_, NA_real_, mean(CF_foldavg_failed, na.rm = TRUE), NA_real_, NA_real_)

        IFE_failure_rate <- c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, mean(IFE_failed, na.rm = TRUE), NA_real_)

        # Summary table

        estimator_names <- c("interaction", "additive", "additive analytical correction", "analytical cf", "analytical cf fold average", "IFE", "TWFE")

        to_save <- rbind(estimator_names, bias, bias_median, cover_rate, sd_est, se_average, cluster_G_average, cluster_C_average, CF_foldavg_failure_rate, IFE_failure_rate)

        rownames(to_save) <- c("Estimator", "Bias", "Median bias", "Coverage", "SD", "Mean SE", "Average G", "Average C", "CF-A failure rate", "IFE failure rate")

        cat("\nN =", N, ", T =", T, ", link =", link, ", theta =", paste0("c(", paste(theta, collapse = ", "), ")"), ", DGP =", DGP, ", nu =", nu, ", rho =", rho, ", kappa =", kappa, "\n")

        print(to_save)

        # Save summary table

        table_name <- paste0("N", N, "analytical", T, "link ", link, "DGP ", DGP, "nu ", nu, "rho ", rho, "kappa ", kappa, "theta ", theta[1], ".csv")

        write.table(to_save, table_name, sep = ";", row.names = TRUE, quote = FALSE)

        # Save all estimates

        estimate_file <- paste0("N", N, "analytical all estimate", T, "link ", link, "DGP ", DGP, "nu", nu, "rho ", rho, "kappa ", kappa, "theta", theta[1], ".csv")

        write.table(estimate, estimate_file, sep = ";", row.names = FALSE, quote = FALSE)

        # Save all estimated SEs

        se_file <- paste0("N", N, "analytical all se", T, "link ", link, "DGP ", DGP, "nu", nu, "rho ", rho, "kappa ", kappa, "theta ", theta[1], ".csv")

        write.table(se, se_file, sep = ";", row.names = FALSE, quote = FALSE)
      }
    }
  }
}
