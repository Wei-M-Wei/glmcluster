# glmcluster
# glmcluster

Replication code for Inference in non-linear panel data after discretizing unobserved heterogeneity

The repository contains the main estimation routines, Monte Carlo simulations, and the empirical application used to compare grouped fixed-effect estimators, cross-fitted estimators, interactive fixed-effect benchmarks, and two-way fixed-effect benchmarks.

## Repository structure

* `functions for nonlinear.R`
  Core estimation and clustering routines. The main functions are `GFE_est()` and `GFE_est_cf()`.

* `simulation replication.R`
  Monte Carlo simulation code for the logit and probit designs.

* `empirical application replication.R`
  Replication code for the empirical application.

* `Sample_use.csv`
  Dataset used by the empirical replication script.

## Main estimators

The replication code reports the following estimators:

* **BC-E**: grouped fixed-effect estimator with analytical bias correction, with the number of clusters selected by the elbow rule.
* **BC-N**: grouped fixed-effect estimator with analytical bias correction, with the number of clusters selected by the variance/noise-floor rule.
* **BC-CFP-E**: pooled cross-fitted bias-corrected estimator using the elbow rule.
* **BC-CFA-E**: fold-averaged cross-fitted bias-corrected estimator using the elbow rule.
* **BC-CFP-N**: pooled cross-fitted bias-corrected estimator using the variance/noise-floor rule.
* **BC-CFA-N**: fold-averaged cross-fitted bias-corrected estimator using the variance/noise-floor rule.
* **IFE-ABC**: interactive fixed-effect estimator with analytical bias correction.
* **TWFE-ABC**: two-way fixed-effect estimator with analytical bias correction.

The grouped estimators allow both additive and interaction specifications and support k-means and k-center clustering.

## Requirements

The scripts use the following R packages:

```r
sandwich
gtools
foreach
parallel
doParallel
doRNG
pcluster
fixest
alpaca
dplyr
NNR (available at https://github.com/Wei-M-Wei/Factor-Bootstrap-replication)
```

Install the required packages before running the replication files.

## Empirical application

The empirical application can be reproduced by running

```r
source("empirical application replication.R")
```

The script expects `Sample_use.csv` to be located in the working directory.

The data are assumed to have the following column structure:

1. unit identifier,
2. time identifier,
3. binary outcome,
4. covariates.

The script checks for missing values, duplicate unit-time observations, and panel balance before estimation. Continuous covariates are standardized, while binary dummy variables are left unscaled.

The empirical replication reports:

1. BC-E,
2. BC-N,
3. BC-CFP-E and BC-CFA-E,
4. BC-CFP-N and BC-CFA-N,
5. IFE-ABC,
6. TWFE-ABC.

## Simulation study

The Monte Carlo experiments can be run with

```r
source("simulation replication.R")
```

The simulation script considers logit and probit models, four DGPs, and several values of \(T\). It reports bias, median bias, empirical coverage, Monte Carlo standard deviation, average estimated standard error, and the average selected numbers of unit and time clusters.

For DGP 4, the parameter is multidimensional and the first coefficient is used as the target parameter.

### Parallel computing

The simulation code is designed for parallel execution. The current script contains

```r
cl <- makeCluster(50)
```

so the number of workers should be adjusted to the available computing resources before running the full Monte Carlo experiment.

### Simulation dependency

The simulation script currently contains

```r
source("DGP truncated.R")
```

Therefore, `DGP truncated.R` must be available in the working directory before running the simulations.

## Core functions

### `GFE_est()`

`GFE_est()` estimates grouped fixed-effect binary-response models.

Important arguments include:

* `link`: binary-response link, such as `"logit"` or `"probit"`;
* `discrete_type`: `"additive"` or `"interaction"`;
* `cluster_type`: `"kmeans"` or `"kcenter"`;
* `group_selection`: `"variance"`, `"elbow"`, `"stability"`, or `"penalized"`;
* `unit_max_groups` and `time_max_groups`: optional upper bounds on the number of clusters.

The function returns the estimated coefficient vector, standard errors, selected unit and time clusters, and, when requested, the analytically bias-corrected estimator.

### `GFE_est_cf()`

`GFE_est_cf()` implements cross-fitting for the grouped fixed-effect estimator.

The routine learns the clustering structure on training folds and assigns the resulting groups to the corresponding evaluation folds. It provides both:

* a **pooled cross-fitted estimator**; and
* a **fold-averaged cross-fitted estimator**.

For the fold-averaged estimator, the reported plug-in standard error is constructed from the fold-specific standard errors.

## Reproducibility

For reproducibility, run the scripts from the repository root so that the source files and data can be found using their relative paths.

The simulation script sets the random seed within each Monte Carlo replication. Results may still depend on the R version, package versions, and numerical optimization routines.

