GFE_est = function(Y, X_list, link, pre_cluster = FALSE, Du_pre, Dv_pre, unit_cluster = NULL, time_cluster = NULL, gamma = 1, discrete_type = 'interaction', cluster_type = 'kmeans', dim_mom_X = NULL, init = 100L, group_selection = 'variance', unit_max_groups = NULL, time_max_groups = NULL, stability_threshold = 0.05, stability_window = 3L, penalty_multiplier = 1, compute_analytical = TRUE, reuse_analytical_fit = FALSE){
  # pre_cluster: if you need to specify the unit and time indicator matrix, set it as 'TRUE'.
  # unit_cluster: if pre_cluster = TRUE, give the specified unit indicator matrix.
  # time_cluster: if pre_cluster = TRUE, give the specified time indicator matrix.
  # gamma: control the choice of the cluster numbers
  # discrete_type: either 'interaction' or 'additive'
  # init: number of clustering initializations
  # group_selection: 'variance', 'elbow', 'stability', or 'penalized'
  # unit_max_groups/time_max_groups: optional upper endpoints, also bounded by
  # floor(2 * N^(1/3)) and floor(2 * T^(1/3)), respectively
  # stability_threshold: largest relative error improvement treated as stable
  # stability_window: number of consecutive future improvements that must be stable
  # penalty_multiplier: multiplier on the BIC-style log(N*T) parameter penalty
  # compute_analytical: set FALSE when only model/se are needed; this skips the
  # additional alpaca fit and analytical bias correction
  # reuse_analytical_fit: for additive models, use the alpaca likelihood fit
  # for both the raw and corrected estimates instead of fitting twice

  N = dim(Y)[1]
  T = dim(Y)[2]
  clusteri <- NULL
  clustert <- NULL
  group_selection <- match.arg(group_selection, c('variance', 'elbow', 'stability', 'penalized'))
  cluster_type <- match.arg(cluster_type, c('kmeans', 'kcenter'))
  automatic_group_selection <- is.null(unit_cluster) && is.null(time_cluster)
  if (!pre_cluster && automatic_group_selection &&
      group_selection == 'elbow' && cluster_type != 'kcenter') {
    stop("group_selection = 'elbow' requires cluster_type = 'kcenter', because the elbow path is defined from the k-center covering radius.")
  }
  if (!pre_cluster && group_selection == 'penalized' && discrete_type != 'additive') {
    stop("group_selection = 'penalized' currently requires discrete_type = 'additive', because its penalty counts the additive fixed-effect parameters.")
  }
  if (reuse_analytical_fit && (discrete_type != 'additive' || !compute_analytical)) {
    stop("reuse_analytical_fit requires discrete_type = 'additive' and compute_analytical = TRUE.")
  }
  if (is.null(dim_mom_X)){
    dim_mom_X = length(X_list)
  }
  if (pre_cluster == FALSE){
    if (cluster_type == 'kmeans'){
      # cluster
      if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
        clusteri <- cluster_general(Y, X_list[seq(dim_mom_X)], N, T, init, type = "long", gamma = gamma, cluster_type = 'kmeans', group_selection = group_selection, max_groups = unit_max_groups, stability_threshold = stability_threshold, stability_window = stability_window, penalty_parameter_increment = T, penalty_sample_size = N * T, penalty_multiplier = penalty_multiplier)
        G <- clusteri$clusters
        klong <- clusteri$res


        clustert <- cluster_general(Y, X_list[seq(dim_mom_X)], N, T, init, type = "tall", gamma = gamma, cluster_type = 'kmeans', group_selection = group_selection, max_groups = time_max_groups, stability_threshold = stability_threshold, stability_window = stability_window, penalty_parameter_increment = N, penalty_sample_size = N * T, penalty_multiplier = penalty_multiplier)
        C <- clustert$clusters
        ktall <- clustert$res
      }else{
        clusteri <- cluster_general(Y, X_list[seq(dim_mom_X)], N, T, init, type = "long" , groups = c(floor(unit_cluster)),  cluster_type = 'kmeans' )
        G <- clusteri$clusters
        klong <- clusteri$res


        clustert <- cluster_general(Y, X_list[seq(dim_mom_X)], N, T, init, type = "tall" , groups = c(floor(time_cluster)),  cluster_type = 'kmeans')
        C <- clustert$clusters
        ktall <- clustert$res
      }
      Du <- matrix(0, N, G)
      Dv <- matrix(0, T, C)

      for (j in seq_len(G)) {
        Du[, j] <- as.numeric(klong$cluster == j)
      }

      for (j in seq_len(C)) {
        Dv[, j] <- as.numeric(ktall$cluster == j)
      }
    }else{
      # cluster
      if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
        clusteri <- cluster_general(Y, X_list[seq(dim_mom_X)], N, T, init, type = "long", gamma = gamma, cluster_type = 'kcenter', group_selection = group_selection, max_groups = unit_max_groups, stability_threshold = stability_threshold, stability_window = stability_window, penalty_parameter_increment = T, penalty_sample_size = N * T, penalty_multiplier = penalty_multiplier)
        G <- clusteri$clusters
        klong <- clusteri$res


        clustert <- cluster_general(Y, X_list[seq(dim_mom_X)], N, T, init, type = "tall", gamma = gamma, cluster_type = 'kcenter', group_selection = group_selection, max_groups = time_max_groups, stability_threshold = stability_threshold, stability_window = stability_window, penalty_parameter_increment = N, penalty_sample_size = N * T, penalty_multiplier = penalty_multiplier)
        C <- clustert$clusters
        ktall <- clustert$res
      }else{
        clusteri <- cluster_general(Y, X_list[seq(dim_mom_X)], N, T, init, type = "long" , groups = c(floor(unit_cluster)),  cluster_type = 'kcenter' )
        G <- clusteri$clusters
        klong <- clusteri$res


        clustert <- cluster_general(Y, X_list[seq(dim_mom_X)], N, T, init, type = "tall" , groups = c(floor(time_cluster)),  cluster_type = 'kcenter')
        C <- clustert$clusters
        ktall <- clustert$res
      }
      Du <- matrix(0, N, G)
      Dv <- matrix(0, T, C)

      for (j in seq_len(G)) {
        Du[, j] <- as.numeric(klong$cluster == j)
      }

      for (j in seq_len(C)) {
        Dv[, j] <- as.numeric(ktall$cluster == j)
      }
    }
  }else if (pre_cluster == TRUE){
    Du = Du_pre
    Dv = Dv_pre
    G = dim(Du)[2]
    C = dim(Dv)[2]
  }

  # Step 1: Expand (i,t) grid with unit i varying fastest ---
  NT <- N * T
  row_index <- rep(1:N, times = T)
  col_index <- rep(1:T, each = N)

  # Step 2: Extract unit and time group labels
  unit_group <- max.col(Du, ties.method = "first")  # length N
  time_group <- max.col(Dv, ties.method = "first")  # length T

  gi <- unit_group[row_index]       # length NT
  ct <- time_group[col_index]      # length NT

  Y_vec <- as.vector(Y)
  X_long <- do.call(cbind, lapply(X_list, as.vector))
  colnames(X_long) <- paste0("X", seq_len(length(X_list)))

  # Step 3: Data frame for feglm
  df <- data.frame(
    Y = Y_vec,
    gi = as.factor(gi),
    ct = as.factor(ct),
    unit = as.factor(row_index),
    time = as.factor(col_index)
  )

  for (k in seq_len(length(X_list))) {
    df[[paste0("X", k)]] <- X_long[, k]
  }
  x_vars <- colnames(X_long)
  fml_interaction <- as.formula(
    paste("Y ~", paste(x_vars, collapse = " + "), "| gi^ct")
  )
  # Step 4: Estimate model with fixed effects
  # Fixed effects: individual x time-cluster (`interaction(unit, ct)`)
  #                group x time (`interaction(gi, time)`)
  if (discrete_type == 'interaction'){
    model <- fixest::feglm(
      fml_interaction,
      data = df,
      family = binomial(link = link),
      fixef.rm = 'none',
      data.save = TRUE,
      glm.iter = 10000,
    )
    se = model$se
    if (compute_analytical) {
      df$gi_ct = interaction(df$gi, df$ct, drop = TRUE)
      fml_interaction_correct <- as.formula(
        paste("Y ~", paste(x_vars, collapse = " + "), "| gi_ct")
      )
      res = alpaca::feglm(fml_interaction_correct,
                          data = df,
                          family = binomial(link = link),
                          control = feglmControl( dev.tol = 1e-8, center.tol = 1e-8,
                                                  iter.max = 10000, drop.pc = FALSE))
      est = biasCorr(res, L = 0, panel.structure = c( "classic"))
      se_ana = summary(est)$cm[,2] # * sqrt((N * T) / ( N*T - N*C - T*G  ))
    } else {
      res <- NULL
      est <- NULL
      se_ana <- NULL
    }
  }else if (discrete_type == 'additive'){
    fml <- as.formula(
      paste("Y ~", paste(x_vars, collapse = " + "), "| unit^ct + gi^time")
    )
    if (!reuse_analytical_fit) {
      model <- fixest::feglm(
        fml,
        data = df,
        family = binomial(link = link),
        fixef.rm = 'none',
        data.save = TRUE,
        glm.iter = 10000
      )
      se = model$se  # * sqrt((N * T) / ( N*T - N*C - T*G  ))
    }
    if (compute_analytical) {
      df$unit_ct <- interaction(df$unit, df$ct, drop = TRUE)
      df$gi_time <- interaction(df$gi,   df$time, drop = TRUE)
      fml <- as.formula(
        paste("Y ~", paste(x_vars, collapse = " + "), "| unit_ct + gi_time")
      )
      res <- alpaca::feglm(
        fml,
        data = df,
        family = binomial(link = link),
        control = feglmControl(
          dev.tol   = 1e-8,
          center.tol= 1e-8,
          iter.max = 10000,
          drop.pc   = FALSE
        )
      )
      est = biasCorr(res, L = 0, panel.structure = c( "classic"))
      se_ana = summary(est)$cm[,2]  # * sqrt((N * T) / ( N*T - N*C - T*G  ))
      if (reuse_analytical_fit) {
        model <- res
        se <- sqrt(diag(solve(res$Hessian)))
      }
    } else {
      res <- NULL
      est <- NULL
      se_ana <- NULL
    }
  }

  return(list(model = model, res_analytical = res, est_analytical = est, se = se, se_ana = se_ana, G = G, C = C, Du = Du, Dv = Dv, unit_group_selection = if (is.null(clusteri)) NULL else clusteri$selection_path, time_group_selection = if (is.null(clustert)) NULL else clustert$selection_path, unit_group_cap = if (pre_cluster || !automatic_group_selection) NULL else clusteri$automatic_max_groups, time_group_cap = if (pre_cluster || !automatic_group_selection) NULL else clustert$automatic_max_groups, unit_theoretical_group_cap = if (pre_cluster || !automatic_group_selection) NULL else clusteri$theoretical_max_groups, time_theoretical_group_cap = if (pre_cluster || !automatic_group_selection) NULL else clustert$theoretical_max_groups, group_selection = if (pre_cluster) 'pre_cluster' else group_selection, stability_threshold = stability_threshold, stability_window = as.integer(stability_window), penalty_multiplier = penalty_multiplier))
}

GFE_est_cf = function(Y, X_list, link, pre_cluster = FALSE, Du_pre, Dv_pre, unit_cluster = NULL, time_cluster = NULL, gamma = 1, discrete_type = 'interaction', cluster_type = 'kmeans', dim_mom_X = NULL, init = 100L, group_selection = 'variance', unit_max_groups = NULL, time_max_groups = NULL, stability_threshold = 0.05, stability_window = 3L, penalty_multiplier = 1, compute_fold_average = TRUE, compute_pooled_uncorrected = TRUE){
  # pre_cluster: if you need to specify the unit and time indicator matrix, set it as 'TRUE'.
  # Du_pre: if pre_cluster = TRUE, give the specified unit indicator list, containing the unit cluster for different folds.
  # Dv_pre: if pre_cluster = TRUE, give the specified time indicator list, containing the time cluster for different folds.
  # unit_cluster: if pre_cluster = FALSE, you can specify the number of unit clusters. unit_cluster is a list containing the unit cluster number of each data-fold.
  # time_cluster: if pre_cluster = FALSE, you can specify the number of time clusters. time_cluster is a list containing the time cluster number of each data-fold.
  # gamma: control the choice of the cluster numbers
  # discrete_type: either 'interaction' or 'additive'
  # group_selection: 'variance', 'elbow', 'stability', or 'penalized'
  # compute_fold_average: set FALSE when only the pooled cross-fitted estimate
  # is needed; this skips the four fold-specific alpaca fits
  # compute_pooled_uncorrected: set FALSE when only the analytically corrected
  # pooled estimate is needed; this skips the additional fixest fit
  # Automatic elbow selection uses k-center covering radii on each training fold.

  N = dim(Y)[1]
  T = dim(Y)[2]
  x_vars <- paste0("X", seq_len(length(X_list)))
  group_selection <- match.arg(group_selection, c('variance', 'elbow', 'stability', 'penalized'))
  cluster_type <- match.arg(cluster_type, c('kmeans', 'kcenter'))
  automatic_group_selection <- is.null(unit_cluster) && is.null(time_cluster)
  if (!pre_cluster && xor(is.null(unit_cluster), is.null(time_cluster))) {
    stop("unit_cluster and time_cluster must either both be NULL or both be supplied.")
  }
  if (!pre_cluster && automatic_group_selection &&
      group_selection == 'elbow' && cluster_type != 'kcenter') {
    stop("group_selection = 'elbow' requires cluster_type = 'kcenter', because each training-fold elbow path is defined from the k-center covering radius.")
  }
  if (!pre_cluster && automatic_group_selection &&
      group_selection == 'penalized' && discrete_type != 'additive') {
    stop("group_selection = 'penalized' currently requires discrete_type = 'additive', because its penalty counts the additive fixed-effect parameters.")
  }
  if (is.null(dim_mom_X)){
    dim_mom_X = length(X_list)
  }
  if (dim_mom_X < 1 || dim_mom_X > length(X_list)){
    stop("dim_mom_X must be between 1 and length(X_list).")
  }
  X_mom_list <- X_list[seq_len(dim_mom_X)]

  # cross-fitting
  folds = 2
  split_indices <- function(dim, folds) {
    indices <- seq_len(dim)
    split(seq_len(dim), cut(indices, folds, labels = FALSE))
  }

  # prepare
  N_folds <- split_indices(N, folds)
  T_folds <- split_indices(T, folds)
  fold_combinations <- expand.grid(N_fold = seq_along(N_folds), T_fold = seq_along(T_folds))
  fold_combinations$fold_id <- (fold_combinations$N_fold - 1) * folds + fold_combinations$T_fold
  fold_combinations <- fold_combinations[order(fold_combinations$fold_id), ]
  rownames(fold_combinations) <- NULL
  fold_count <- nrow(fold_combinations)
  beta_estimate_cf = rep(0, nrow(fold_combinations))
  se_cf = rep(0, nrow(fold_combinations))
  df_list = as.list(rep(0, fold_count))
  N_list = as.list(rep(0, fold_count))
  T_list = as.list(rep(0, fold_count))
  eta_cf = as.list(rep(0, fold_count))
  unit_cluster_all = as.list(rep(0, fold_count))
  time_cluster_all = as.list(rep(0, fold_count))
  unit_group_selection_all <- vector("list", fold_count)
  time_group_selection_all <- vector("list", fold_count)
  unit_group_cap_by_fold <- rep(NA_integer_, fold_count)
  time_group_cap_by_fold <- rep(NA_integer_, fold_count)
  unit_theoretical_group_cap_by_fold <- rep(NA_integer_, fold_count)
  time_theoretical_group_cap_by_fold <- rep(NA_integer_, fold_count)
  G_by_fold <- integer(fold_count)
  C_by_fold <- integer(fold_count)
  training_fold_by_evaluation <- rep(NA_integer_, fold_count)
  beta_fold = matrix(NA_real_, fold_count, length(x_vars))
  colnames(beta_fold) <- x_vars
  se_fold = matrix(NA_real_, fold_count, length(x_vars))
  colnames(se_fold) <- x_vars
  freedom = 0
  G_av = 0
  C_av = 0

  fold_pair <- setNames(fold_count + 1 - seq_len(fold_count), seq_len(fold_count))

  resolve_fold_groups <- function(value, fold_index, argument_name) {
    if (is.list(value)) {
      if (length(value) == 1L) {
        value <- value[[1L]]
      } else if (length(value) == fold_count) {
        value <- value[[fold_index]]
      } else {
        stop(argument_name, " must have length one or one entry per evaluation fold.")
      }
    } else if (length(value) == fold_count) {
      value <- value[fold_index]
    }
    if (length(value) != 1L || !is.finite(value) || value < 1) {
      stop(argument_name, " must provide one finite positive group count per fold.")
    }
    as.integer(floor(value))
  }

  compute_long_training <- function(Y_train, X_train_list, cc = 0) {
    N_train <- nrow(Y_train)
    T_train <- ncol(Y_train)
    dimtheta <- length(X_train_list)
    mdim <- dimtheta + 1
    mom_i <- rowMeans(Y_train)
    for (X in X_train_list) mom_i <- cbind(mom_i, rowMeans(X))
    mom_micro <- Y_train
    for (X in X_train_list) mom_micro <- cbind(mom_micro, X)
    col_mean <- numeric(mdim)
    col_sd <- numeric(mdim)
    for (j in seq_len(mdim)) {
      col_mean[j] <- mean(mom_i[, j])
      col_sd[j] <- sd(mom_i[, j])
      if (col_sd[j] > 0) {
        cols <- ((j - 1) * T_train + 1):(j * T_train)
        mom_micro[, cols] <- (mom_micro[, cols] - col_mean[j]) / col_sd[j]
        mom_i[, j] <- (mom_i[, j] - col_mean[j]) / col_sd[j]
      }
    }
    if (cc == 1) {
      var_w <- var_tot <- matrix(0, mdim, mdim)
      for (j in seq_len(mdim)) {
        for (k in seq_len(mdim)) {
          cols_j <- ((j - 1) * T_train + 1):(j * T_train)
          cols_k <- ((k - 1) * T_train + 1):(k * T_train)
          var_w[j, k] <- mean(colMeans((mom_micro[, cols_j] - mom_i[, j]) *
                                         (mom_micro[, cols_k] - mom_i[, k]))) / T_train
          var_tot[j, k] <- mean(mom_i[, j] * mom_i[, k])
        }
      }
      var_b <- var_tot - var_w
      mat_scali <- diag(diag(var_b)) / diag(diag(var_tot))
      mat_scali <- diag(diag(mat_scali))
      mat_scali <- pmax(mat_scali, 0)
      if (any(is.nan(mat_scali)) || any(is.infinite(mat_scali)) || max(mat_scali) == 0) {
        mat_scali <- diag(mdim)
      }
    } else {
      mat_scali <- diag(mdim)
    }
    moment_weights <- diag(mat_scali)
    mom_i <- sweep(mom_i, 2, moment_weights, `*`)
    mom_micro <- sweep(mom_micro, 2, rep(moment_weights, each = T_train), `*`)
    dimnames(mom_i) <- NULL
    dimnames(mom_micro) <- NULL
    expanded_mom_i <- mom_i[, rep(seq_len(mdim), each = T_train), drop = FALSE]
    variance <- sum((mom_micro - expanded_mom_i)^2) / (N_train * T_train^2)
    list(data = mom_i, variance = variance, dim_size = N_train,
         col_mean = col_mean, col_sd = col_sd, mat_scali = mat_scali)
  }

  compute_long_apply <- function(Y_target, X_target_list, col_mean, col_sd, mat_scali) {
    T_target <- ncol(Y_target)
    dimtheta <- length(X_target_list)
    mdim <- dimtheta + 1
    mom_i <- rowMeans(Y_target)
    for (X in X_target_list) mom_i <- cbind(mom_i, rowMeans(X))
    mom_micro <- Y_target
    for (X in X_target_list) mom_micro <- cbind(mom_micro, X)
    for (j in seq_len(mdim)) {
      if (col_sd[j] > 0) {
        cols <- ((j - 1) * T_target + 1):(j * T_target)
        mom_micro[, cols] <- (mom_micro[, cols] - col_mean[j]) / col_sd[j]
        mom_i[, j] <- (mom_i[, j] - col_mean[j]) / col_sd[j]
      }
    }
    mom_i <- mom_i %*% mat_scali
    mom_i
  }

  compute_tall_training <- function(Y_train, X_train_list, cc = 0) {
    N_train <- nrow(Y_train)
    T_train <- ncol(Y_train)
    dimtheta <- length(X_train_list)
    mdim <- dimtheta + 1
    mom_t <- colMeans(Y_train)
    for (X in X_train_list) mom_t <- cbind(mom_t, colMeans(X))
    mom_micro2 <- t(Y_train)
    for (X in X_train_list) mom_micro2 <- cbind(mom_micro2, t(X))
    col_mean <- mean(mom_t[, 1])
    col_sd <- sd(mom_t[, 1])
    if (col_sd > 0) {
      cols <- 1:N_train
      mom_micro2[, cols] <- (mom_micro2[, cols] - col_mean) / col_sd
      mom_t[, 1] <- (mom_t[, 1] - col_mean) / col_sd
    }
    if (cc == 1) {
      var_w2 <- var_tot2 <- matrix(0, mdim, mdim)
      for (j in seq_len(mdim)) {
        for (k in seq_len(mdim)) {
          cols_j <- ((j - 1) * N_train + 1):(j * N_train)
          cols_k <- ((k - 1) * N_train + 1):(k * N_train)
          var_w2[j, k] <- mean((mom_micro2[, cols_j] - mom_t[, j]) *
                                 (mom_micro2[, cols_k] - mom_t[, k])) / N_train
          var_tot2[j, k] <- mean(mom_t[, j] * mom_t[, k])
        }
      }
      var_b2 <- var_tot2 - var_w2
      mat_scali2 <- diag(diag(var_b2)) / diag(diag(var_tot2))
      mat_scali2 <- diag(diag(mat_scali2))
      mat_scali2 <- pmax(mat_scali2, 0)
    } else {
      mat_scali2 <- diag(mdim)
    }
    if (any(is.nan(mat_scali2)) || any(is.infinite(mat_scali2)) || max(mat_scali2) == 0) {
      mat_scali2 <- diag(mdim)
    }
    moment_weights <- diag(mat_scali2)
    mom_t <- sweep(mom_t, 2, moment_weights, `*`)
    mom_micro2 <- sweep(mom_micro2, 2, rep(moment_weights, each = N_train), `*`)
    dimnames(mom_t) <- NULL
    dimnames(mom_micro2) <- NULL
    expanded_mom_t <- mom_t[, rep(seq_len(mdim), each = N_train), drop = FALSE]
    variance <- sum((mom_micro2 - expanded_mom_t)^2) / (T_train * N_train^2)
    list(data = mom_t, variance = variance, dim_size = T_train,
         col_mean = col_mean, col_sd = col_sd, mat_scali = mat_scali2)
  }

  compute_tall_apply <- function(Y_target, X_target_list, col_mean, col_sd, mat_scali) {
    N_target <- nrow(Y_target)
    dimtheta <- length(X_target_list)
    mdim <- dimtheta + 1
    mom_t <- colMeans(Y_target)
    for (X in X_target_list) mom_t <- cbind(mom_t, colMeans(X))
    mom_micro2 <- t(Y_target)
    for (X in X_target_list) mom_micro2 <- cbind(mom_micro2, t(X))
    if (col_sd > 0) {
      cols <- 1:N_target
      mom_micro2[, cols] <- (mom_micro2[, cols] - col_mean) / col_sd
      mom_t[, 1] <- (mom_t[, 1] - col_mean) / col_sd
    }
    mom_t <- mom_t %*% mat_scali
    mom_t
  }

  assign_to_centers <- function(data, centers) {
    data <- as.matrix(data)
    centers <- as.matrix(centers)
    if (nrow(centers) == 1) return(rep(1, nrow(data)))
    dists <- sapply(seq_len(nrow(centers)), function(k) {
      rowSums((data - matrix(centers[k, ], nrow(data), ncol(data), byrow = TRUE))^2)
    })
    max.col(-dists, ties.method = "first")
  }

  for (i in seq_len(nrow(fold_combinations))) {

    N_idx <- N_folds[[fold_combinations$N_fold[i]]]
    T_idx <- T_folds[[fold_combinations$T_fold[i]]]

    N_list[[i]] = N_idx
    T_list[[i]] = T_idx

    dim_N <- length(N_idx)
    dim_T <- length(T_idx)

    Y_comp_long <- Y[N_idx, -T_idx]
    X_comp_long_list <- lapply(X_list, function(X) X[N_idx, -T_idx])

    Y_comp_tall <- Y[-N_idx, T_idx]
    X_comp_tall_list <- lapply(X_list, function(X) X[-N_idx, T_idx])

    long_fit <- NULL
    tall_fit <- NULL
    train_row <- NA_integer_
    if (pre_cluster == FALSE){
      train_row <- fold_pair[[as.character(i)]]
      training_fold_by_evaluation[i] <- train_row
      N_train <- N_folds[[fold_combinations$N_fold[train_row]]]
      T_train <- T_folds[[fold_combinations$T_fold[train_row]]]
      Y_train <- Y[N_train, T_train, drop = FALSE]
      X_train_mom_list <- lapply(X_mom_list, function(X) X[N_train, T_train, drop = FALSE])
      Y_target <- Y[N_idx, T_idx, drop = FALSE]
      X_target_mom_list <- lapply(X_mom_list, function(X) X[N_idx, T_idx, drop = FALSE])
      unit_groups_i <- if (automatic_group_selection) NULL else {
        resolve_fold_groups(unit_cluster, i, "unit_cluster")
      }
      time_groups_i <- if (automatic_group_selection) NULL else {
        resolve_fold_groups(time_cluster, i, "time_cluster")
      }

      long_fit <- cluster_general(
        Y_train,
        X_train_mom_list,
        length(N_train),
        length(T_train),
        init,
        type = "long",
        groups = unit_groups_i,
        gamma = gamma,
        cluster_type = cluster_type,
        group_selection = group_selection,
        max_groups = unit_max_groups,
        stability_threshold = stability_threshold,
        stability_window = stability_window,
        penalty_parameter_increment = length(T_train),
        penalty_sample_size = length(N_train) * length(T_train),
        penalty_multiplier = penalty_multiplier
      )
      G <- long_fit$clusters
      target_long <- compute_long_training(Y_target, X_target_mom_list)$data
      du_assign <- assign_to_centers(target_long, long_fit$centers)
      Du <- matrix(0, dim_N, G)
      for (j in seq_len(G)) {
        Du[, j] <- as.numeric(du_assign == j)
      }

      tall_fit <- cluster_general(
        Y_train,
        X_train_mom_list,
        length(N_train),
        length(T_train),
        init,
        type = "tall",
        groups = time_groups_i,
        gamma = gamma,
        cluster_type = cluster_type,
        group_selection = group_selection,
        max_groups = time_max_groups,
        stability_threshold = stability_threshold,
        stability_window = stability_window,
        penalty_parameter_increment = length(N_train),
        penalty_sample_size = length(N_train) * length(T_train),
        penalty_multiplier = penalty_multiplier
      )
      C <- tall_fit$clusters
      target_tall <- compute_tall_training(Y_target, X_target_mom_list)$data
      dv_assign <- assign_to_centers(target_tall, tall_fit$centers)
      Dv <- matrix(0, dim_T, C)
      for (j in seq_len(C)) {
        Dv[, j] <- as.numeric(dv_assign == j)
      }
      if (automatic_group_selection) {
        unit_group_cap_by_fold[i] <- long_fit$automatic_max_groups
        time_group_cap_by_fold[i] <- tall_fit$automatic_max_groups
        unit_theoretical_group_cap_by_fold[i] <- long_fit$theoretical_max_groups
        time_theoretical_group_cap_by_fold[i] <- tall_fit$theoretical_max_groups
      }
    }else if (pre_cluster == TRUE){
      Du = Du_pre[[i]]
      Dv = Dv_pre[[i]]
      G = dim(Du)[2]
      C = dim(Dv)[2]
    }

    G_by_fold[i] <- G
    C_by_fold[i] <- C
    if (!is.null(long_fit) && !is.null(long_fit$selection_path)) {
      unit_group_selection_all[[i]] <- long_fit$selection_path
      unit_group_selection_all[[i]]$evaluation_fold <- i
      unit_group_selection_all[[i]]$training_fold <- train_row
    }
    if (!is.null(tall_fit) && !is.null(tall_fit$selection_path)) {
      time_group_selection_all[[i]] <- tall_fit$selection_path
      time_group_selection_all[[i]]$evaluation_fold <- i
      time_group_selection_all[[i]]$training_fold <- train_row
    }

    unit_cluster_all[[i]] = Du
    time_cluster_all[[i]] = Dv


    # Step 1: Expand (i,t) grid with unit i varying fastest ---
    NT_idx <- dim_N * dim_T
    row_index <- rep(seq_len(dim_N), times = dim_T)
    col_index <- rep(seq_len(dim_T), each = dim_N)
    unit_index <- rep(N_idx, times = dim_T)
    time_index <- rep(T_idx, each = dim_N)

    # Step 2: Extract unit and time group labels
    unit_group <- max.col(Du, ties.method = "first")  # length N
    time_group <- max.col(Dv, ties.method = "first")  # length T

    gi <- unit_group[row_index]       # length NT
    ct <- time_group[col_index]      # length NT

    Y_vec <- as.vector(Y[N_idx,T_idx])
    X_long <- do.call(
      cbind,
      lapply(seq_along(X_list), function(i) {
        as.vector(X_list[[i]][N_idx, T_idx])
      })
    )
    colnames(X_long) <- paste0("X", seq_len(length(X_list)))

    # Step 3: Data frame for feglm
    df <- data.frame(
      Y = Y_vec,
      gi = as.factor(gi),
      ct = as.factor(ct),
      unit = as.factor(unit_index),
      time = as.factor(time_index)
    )

    for (k in seq_len(length(X_list))) {
      df[[paste0("X", k)]] <- X_long[, k]
    }
    if (compute_fold_average && discrete_type == 'interaction'){
      df$gi_ct <- interaction(df$gi, df$ct, drop = TRUE)
      fml_fold_correct <- as.formula(
        paste("Y ~", paste(x_vars, collapse = " + "), "| gi_ct")
      )
      res_fold <- alpaca::feglm(
        fml_fold_correct,
        data = df,
        family = binomial(link = link),
        control = feglmControl(
          dev.tol   = 1e-6,
          center.tol= 1e-6,
          iter.max = 1000,
          drop.pc   = FALSE
        )
      )
      est_fold <- biasCorr(res_fold, L = 0, panel.structure = c("classic"))
    }else if (compute_fold_average && discrete_type == 'additive'){
      df$unit_ct <- interaction(df$unit, df$ct, drop = TRUE)
      df$gi_time <- interaction(df$gi, df$time, drop = TRUE)
      fml_fold_correct <- as.formula(
        paste("Y ~", paste(x_vars, collapse = " + "), "| unit_ct + gi_time")
      )
      res_fold <- alpaca::feglm(
        fml_fold_correct,
        data = df,
        family = binomial(link = link),
        control = feglmControl(
          dev.tol   = 1e-6,
          center.tol= 1e-6,
          iter.max = 1000,
          drop.pc   = FALSE
        )
      )
      est_fold <- biasCorr(res_fold, L = 0, panel.structure = c("classic"))
    }
    if (compute_fold_average) {
      beta_i <- coef(est_fold)
      idx <- intersect(names(beta_i), x_vars)
      if (length(idx) > 0) {
        beta_fold[i, idx] <- beta_i[idx]
      }
      se_cm <- summary(est_fold)$cm
      if (is.matrix(se_cm)) {
        se_vals <- se_cm[, 2]
        names(se_vals) <- rownames(se_cm)
      } else {
        se_vals <- se_cm[2]
      }
      if (is.null(names(se_vals)) && !is.null(names(beta_i))) {
        names(se_vals) <- names(beta_i)
      }
      if (length(idx) > 0) {
        se_fold[i, idx] <- se_vals[idx]
      }
    }
    df_fold = df
    df_fold$fold <- i
    df_list[[i]] <- df_fold

    freedom = freedom + dim_N * dim_T - dim_T * G - dim_N * C
    G_av = G_av + G
    C_av = C_av + C
  }

  # another way to get the cross-fitting eatimate
  df_all <- do.call(rbind, df_list)
  for(d in seq_len(fold_count)){
    df_all[[paste0("ct_", d)]] <- ifelse(df_all$fold == d,
                                         as.character(df_all$ct),
                                         "0")
    df_all[[paste0("gi_", d)]] <- ifelse(df_all$fold == d,
                                         as.character(df_all$gi),
                                         "0")
  }
  # make sure they are factors
  df_all[ , grep("^(ct_|gi_)", names(df_all))] <- lapply(
    df_all[ , grep("^(ct_|gi_)", names(df_all))], factor
  )

  if (compute_fold_average) {
    beta_fold_avg <- colMeans(beta_fold, na.rm = TRUE)
    n_valid_folds <- colSums(!is.na(se_fold))
    se_fold_avg <- sqrt(colSums(se_fold^2, na.rm = TRUE)) / n_valid_folds
    beta_fold_sd <- apply(beta_fold, 2, sd, na.rm = TRUE)
  } else {
    beta_fold_avg <- setNames(rep(NA_real_, length(x_vars)), x_vars)
    se_fold_avg <- setNames(rep(NA_real_, length(x_vars)), x_vars)
    beta_fold_sd <- setNames(rep(NA_real_, length(x_vars)), x_vars)
  }
  if (discrete_type == 'interaction'){
    fold_terms <- paste0("gi_", seq_len(fold_count), "^ct_", seq_len(fold_count))
    fml_cf <- as.formula(
      paste("Y ~", paste(x_vars, collapse = " + "), "|", paste(fold_terms, collapse = " + "))
    )
    if (compute_pooled_uncorrected) {
      model = fixest::feglm(
        fml_cf,
        data = df_all,
        family = binomial(link = link),
        fixef.rm = "none",
        data.save = TRUE,
        glm.iter = 1000
      )
      se = model$se
    } else {
      model <- NULL
      se <- NULL
    }
    df_all$gi_ct_fold <- interaction(df_all$gi, df_all$ct, df_all$fold, drop = TRUE)
    fml_cf_correct <- as.formula(
      paste("Y ~", paste(x_vars, collapse = " + "), "| gi_ct_fold")
    )
    res <- alpaca::feglm(
      fml_cf_correct,
      data = df_all,
      family = binomial(link = link),
      control = feglmControl(
        dev.tol   = 1e-6,
        center.tol= 1e-6,
        iter.max = 1000,
        drop.pc   = FALSE
      )
    )
    est = biasCorr(res, L = 0, panel.structure = c( "classic"))
    se_ana = summary(est)$cm[,2]
  }else if (discrete_type == 'additive'){
    fml_cf <- as.formula(
      paste("Y ~", paste(x_vars, collapse = " + "), "| interaction(unit, ct, fold) + interaction(gi, time, fold)")
    )
    if (compute_pooled_uncorrected) {
      model = fixest::feglm(
        fml_cf,
        data = df_all,
        family = binomial(link = link),
        fixef.rm = "none",
        data.save = TRUE,
        glm.tol = 1e-6,
        glm.iter = 1000
      )
      se = model$se
    } else {
      model <- NULL
      se <- NULL
    }
    df_all$unit_ct_fold <- interaction(df_all$unit, df_all$ct, df_all$fold, drop = TRUE)
    df_all$gi_time_fold <- interaction(df_all$gi,   df_all$time, df_all$fold, drop = TRUE)
    fml_cf_correct <- as.formula(
      paste("Y ~", paste(x_vars, collapse = " + "), "| unit_ct_fold + gi_time_fold")
    )
    res <- alpaca::feglm(
      fml_cf_correct,
      data = df_all,
      family = binomial(link = link),
      control = feglmControl(
        dev.tol   = 1e-6,
        center.tol= 1e-6,
        iter.max = 1000,
        drop.pc   = FALSE
      )
    )
    est = biasCorr(res, L = 0, panel.structure = c( "classic"))
    se_ana = summary(est)$cm[,2]
  }
  return(list(model = model, res_analytical = res, est_analytical = est, se = se, se_ana = se_ana,
              beta_fold = beta_fold, se_fold = se_fold, beta_fold_avg = beta_fold_avg, se_fold_avg = se_fold_avg,
              beta_fold_sd = beta_fold_sd,
              G = G_av/fold_count, C = C_av/fold_count, unit_cluster_all = unit_cluster_all, time_cluster_all = time_cluster_all,
              G_by_fold = G_by_fold, C_by_fold = C_by_fold,
              unit_group_selection_all = unit_group_selection_all,
              time_group_selection_all = time_group_selection_all,
              unit_group_cap_by_fold = unit_group_cap_by_fold,
              time_group_cap_by_fold = time_group_cap_by_fold,
              unit_theoretical_group_cap_by_fold = unit_theoretical_group_cap_by_fold,
              time_theoretical_group_cap_by_fold = time_theoretical_group_cap_by_fold,
              training_fold_by_evaluation = training_fold_by_evaluation,
              group_selection = if (pre_cluster) 'pre_cluster' else if (automatic_group_selection) group_selection else 'fixed',
              stability_threshold = stability_threshold,
              stability_window = as.integer(stability_window),
              penalty_multiplier = penalty_multiplier,
              N_list = N_list, T_list = T_list))
}
# GFE_neyman = function(Y, X_list, X_list_weighted, link, weights = rep(1,length(Y)), pre_cluster = FALSE, Du_pre, Dv_pre, unit_cluster = NULL, time_cluster = NULL, gamma = 1, discrete_type = 'interaction'){
#
#   N = dim(Y)[1]
#   T = dim(Y)[2]
#
#   if (pre_cluster == FALSE){
#     # cluster
#     if (is.null(unit_cluster) == 1 & is.null(time_cluster) == 1){
#       clusteri <- cluster_general(Y, X_list_weighted, N, T, init, type = "long", gamma = gamma  )
#       G <- clusteri$clusters
#       klong <- clusteri$res
#
#
#       clustert <- cluster_general(Y, X_list_weighted, N, T, init, type = "tall" , gamma = gamma )
#       C <- clustert$clusters
#       ktall <- clustert$res
#     }else{
#       clusteri <- cluster_general(Y, X_list_weighted, N, T, init, type = "long" , groups = c(floor(unit_cluster)) )
#       G <- clusteri$clusters
#       klong <- clusteri$res
#
#
#       clustert <- cluster_general(Y, X_list_weighted, N, T, init, type = "tall" , groups = c(floor(time_cluster)))
#       C <- clustert$clusters
#       ktall <- clustert$res
#     }
#     Du <- matrix(0, N, G)
#     Dv <- matrix(0, T, C)
#
#     for (j in seq_len(G)) {
#       Du[, j] <- as.numeric(klong$cluster == j)
#     }
#
#     for (j in seq_len(C)) {
#       Dv[, j] <- as.numeric(ktall$cluster == j)
#     }
#   }else{
#     Du = Du_pre
#     Dv = Dv_pre
#   }
#
#   # Step 1: Expand (i,t) grid with unit i varying fastest ---
#   NT <- N * T
#   row_index <- rep(1:N, times = T)
#   col_index <- rep(1:T, each = N)
#
#   # Step 2: Extract unit and time group labels
#   unit_group <- apply(Du, 1, function(row) which(row == 1))  # length N
#   time_group <- apply(Dv, 1, function(row) which(row == 1))  # length T
#
#   gi <- unit_group[row_index]       # length NT
#   ct <- time_group[col_index]      # length NT
#
#   Y_vec <- as.vector(Y)
#   X_long <- as.matrix(unlist(lapply(X_list, as.vector)))
#
#   # Step 3: Data frame for feglm
#   df <- data.frame(
#     Y = Y_vec,
#     gi = as.factor(gi),
#     ct = as.factor(ct),
#     unit = as.factor(row_index),
#     time = as.factor(col_index)
#   )
#
#   for (k in 1:length(X_list)) {
#     df[[paste0("X", k)]] <- X_long[, k]
#   }
#
#   # Step 4: Estimate probit model with fixed effects
#   # Fixed effects: individual x time-cluster (`interaction(unit, ct)`)
#   #                group x time (`interaction(gi, time)`)
#
#   if (discrete_type == 'interaction'){
#     model <- fixest::feglm(
#       Y ~ X1 | interaction(gi, ct),
#       data = df,
#       weights= weights,
#       family = binomial(link = link),
#       fixef.rm = 'none',
#       data.save = TRUE,
#       glm.iter = 200
#     )
#   }else if (discrete_type == 'additive'){
#     model <- fixest::feglm(
#       Y ~ X1 | interaction(unit, ct) + interaction(gi, time),
#       data = df,
#       weights= weights,
#       family = binomial(link = link),
#       fixef.rm = 'none',
#       data.save = TRUE,
#       glm.iter = 200
#     )
#   }
#
#   return(list(model = model, G = G, C = C, Du = Du, Dv = Dv))
# }

# GFE_neyman_est = function(Y, X_list, link, pre_cluster = FALSE, Du_pre, Dv_pre, unit_cluster = NULL, time_cluster = NULL, gamma = 1, discrete_type = 'interaction'){
#
#   N = dim(Y)[1]
#   T = dim(Y)[2]
#
#   model = GFE_est(Y = Y, X_list = X_list, link = link, pre_cluster = pre_cluster, Du_pre = Du_pre, Dv_pre = Dv_pre, unit_cluster = unit_cluster, time_cluster = time_cluster, gamma = gamma, discrete_type = discrete_type )
#   G = model$G
#   C = model$C
#   Du = model$Du
#   Dv = model$Dv
#   eta <- predict(model$model, type = "link")
#
#   if (model$model$family$link  =='logit'){
#     G_hat = plogis(eta)
#     w = dlogis(eta)
#     sigma = sqrt(G_hat*(1-G_hat))
#     f = w/sigma
#     f_matrix = list(matrix(f, N,T))
#   }else{
#     G_hat = pnorm(eta)
#     w = dnorm(eta)
#     sigma = sqrt(G_hat*(1-G_hat))
#     f = w/sigma
#     f_matrix = list(matrix(f, N,T))
#   }
#
#   weighted_X = Map(`*`, list(matrix(w, N, T)), X_list)
#   model_weighted <- GFE_neyman(Y = Y, X_list = X_list, X_list_weighted = X_list, link = link, pre_cluster = pre_cluster, Du_pre = Du_pre, Dv_pre = Dv_pre, unit_cluster = unit_cluster, time_cluster = time_cluster, gamma = gamma, discrete_type = discrete_type )
#
#
#   return(list(model = model_weighted$model, G = model_weighted$G, C = model_weighted$C, Du = model_weighted$Du, Dv = model_weighted$Dv))
# }

# cluster_general <- function(Y, X_list, N, T, init = 30, type = "long", groups = NULL) {
#   if (type == "long") {
#     Y_mean <- rowMeans(Y)
#     X_means <- lapply(X_list, rowMeans)
#     data <- cbind(Y_mean, do.call(cbind, X_means))
#     variance <- sum(sapply(1:T, function(t) {
#       norm(cbind(Y[, t], sapply(X_list, function(X) X[, t])) - data, type = "F")^2
#     })) / (N * T^2)
#     dim_size <- N
#   } else if (type == "tall") {
#     Y_mean <- colMeans(Y)
#     X_means <- lapply(X_list, colMeans)
#     data <- cbind(Y_mean, do.call(cbind, X_means))
#     variance <- sum(sapply(1:N, function(i) {
#       norm(cbind(Y[i, ], sapply(X_list, function(X) X[i, ])) - data, type = "F")^2
#     })) / (N^2 * T)
#     dim_size <- T
#   } else {
#     stop("Invalid type. Use 'long' for rows or 'tall' for columns.")
#   }
#
#   # Clustering logic
#   if (!is.null(groups)) {
#     clusters <- groups
#     k_result <- kmeans(data, centers = clusters, algorithm = "Hartigan-Wong", nstart = init)
#   } else {
#     clusters <- 1
#     repeat {
#       k_result <- kmeans(data, centers = clusters, algorithm = "Hartigan-Wong", nstart = init)
#       if (k_result$tot.withinss / dim_size <= variance) break
#       clusters <- clusters + 1
#     }
#   }
#   list(res = k_result, clusters = clusters)
# }

reshape_to_matrix <- function(data, id_var, time_var, value_var) {
  reshaped <- reshape(data[, c(id_var, time_var, value_var)],
                      idvar = id_var,
                      timevar = time_var,
                      direction = "wide")
  as.matrix(reshaped[, -1])  # Drop the 'id' column and return the matrix
}

split_panel = function(Y, X_list, model, link, weights = rep(1,length(Y))){
  index_1 = split_groups_balanced(model$Du)$indices[[1]]
  index_2 = split_groups_balanced(model$Du)$indices[[2]]

  X_long_1 <- do.call(cbind, lapply(X_list, function(X) as.vector(X[index_1, ])))
  colnames(X_long_1) <- paste0("X", seq_len(length(X_list)))
  df1 = data.frame(
    Y = as.vector(Y[index_1,]),
    unit = as.factor(index_1),
    time = as.factor(rep(1:T, each = length(index_1)))
  )
  for (k in seq_len(length(X_list))) {
    df1[[paste0("X", k)]] <- X_long_1[, k]
  }

  X_long_2 <- do.call(cbind, lapply(X_list, function(X) as.vector(X[index_2, ])))
  colnames(X_long_2) <- paste0("X", seq_len(length(X_list)))
  df2 = data.frame(
    Y = as.vector(Y[index_2,]),
    unit = as.factor(index_2),
    time = as.factor(rep(1:T, each = length(index_2)))
  )
  for (k in seq_len(length(X_list))) {
    df2[[paste0("X", k)]] <- X_long_2[, k]
  }
  x_vars <- paste0("X", seq_len(length(X_list)))
  fml_split <- as.formula(paste("Y ~", paste(x_vars, collapse = " + "), "| unit + time"))

  estimator_1 = fixest::feglm(
    fml_split,
    data = df1,
    family = binomial(link = link),
    fixef.rm = 'none',
    data.save = TRUE
  )
  estimator_2 = fixest::feglm(
    fml_split,
    data = df2,
    family = binomial(link = link),
    fixef.rm = 'none',
    data.save = TRUE
  )
  estimator_unit = (coef(estimator_1) + coef(estimator_2))/2

  index_1 = split_groups_balanced(model$Dv)$indices[[1]]
  index_2 = split_groups_balanced(model$Dv)$indices[[2]]
  X_long_1 <- do.call(cbind, lapply(X_list, function(X) as.vector(X[, index_1])))
  colnames(X_long_1) <- paste0("X", seq_len(length(X_list)))
  df1 = data.frame(
    Y = as.vector(Y[,index_1]),
    unit = as.factor(rep(1:N, each = length(index_1))),
    time = as.factor(index_1)
  )
  for (k in seq_len(length(X_list))) {
    df1[[paste0("X", k)]] <- X_long_1[, k]
  }

  X_long_2 <- do.call(cbind, lapply(X_list, function(X) as.vector(X[, index_2])))
  colnames(X_long_2) <- paste0("X", seq_len(length(X_list)))
  df2 = data.frame(
    Y = as.vector(Y[,index_2]),
    unit = as.factor(rep(1:N, each = length(index_2))),
    time = as.factor(index_2)
  )
  for (k in seq_len(length(X_list))) {
    df2[[paste0("X", k)]] <- X_long_2[, k]
  }
  if (length(unique(df1$Y)) == 1){
    estimator_time_1 = 0
  }else{
    estimator_time_1 = fixest::feglm(
      fml_split,
      data = df1,
      family = binomial(link = link),
      fixef.rm = 'none',
      data.save = TRUE
    )
    estimator_time_1= coef(estimator_time_1)
  }
  if (length(unique(df2$Y)) == 1){
    estimator_time_2 = 0
  }else{
    estimator_time_2 = fixest::feglm(
      fml_split,
      data = df2,
      family = binomial(link = link),
      fixef.rm = 'none',
      data.save = TRUE
    )
    estimator_time_2= coef(estimator_time_2)
  }
  estimator_time = (estimator_time_1 + estimator_time_2)/2

  return(list(estimator_time = estimator_time, estimator_unit = estimator_unit))

}

select_elbow_groups <- function(group_counts, clustering_error) {
  group_counts <- as.integer(group_counts)
  clustering_error <- as.numeric(clustering_error)

  if (length(group_counts) != length(clustering_error) || length(group_counts) == 0) {
    stop("group_counts and clustering_error must have the same positive length.")
  }
  if (any(!is.finite(clustering_error)) || any(clustering_error < 0)) {
    stop("clustering_error must contain finite non-negative values.")
  }
  if (is.unsorted(group_counts, strictly = TRUE)) {
    stop("group_counts must be strictly increasing.")
  }

  monotone_error <- cummin(clustering_error)
  error_drop <- monotone_error[1] - monotone_error[length(monotone_error)]
  error_scale <- max(1, abs(monotone_error[1]))

  normalized_groups <- if (length(group_counts) == 1) {
    0
  } else {
    (group_counts - group_counts[1]) /
      (group_counts[length(group_counts)] - group_counts[1])
  }
  normalized_error <- rep(NA_real_, length(group_counts))
  flat_path <- error_drop <= sqrt(.Machine$double.eps) * error_scale

  if (length(group_counts) == 1 || flat_path) {
    elbow_score <- rep(0, length(group_counts))
    selected_index <- 1L
  } else {
    normalized_error <- (monotone_error - monotone_error[length(monotone_error)]) /
      error_drop
    elbow_score <- (1 - normalized_groups) - normalized_error
    elbow_score[c(1, length(elbow_score))] <- 0
    # which.max() returns the first maximizer, implementing the parsimonious
    # tie rule in favour of the smallest candidate group count.
    selected_index <- which.max(elbow_score)
  }

  relative_improvement <- rep(NA_real_, length(monotone_error))
  if (length(monotone_error) > 1) {
    previous_error <- monotone_error[-length(monotone_error)]
    relative_improvement[-1] <- ifelse(
      previous_error > 0,
      (previous_error - monotone_error[-1]) / previous_error,
      0
    )
  }

  path <- data.frame(
    groups = group_counts,
    clustering_error = clustering_error,
    monotone_error = monotone_error,
    normalized_group_count = normalized_groups,
    normalized_clustering_error = normalized_error,
    relative_improvement = relative_improvement,
    elbow_score = elbow_score,
    flat_path = flat_path,
    selected = seq_along(group_counts) == selected_index,
    stringsAsFactors = FALSE
  )

  list(groups = group_counts[selected_index], path = path)
}

# Select the first group count followed by a sustained plateau in the error path.
select_stability_groups <- function(group_counts, clustering_error,
                                    stability_threshold = 0.05,
                                    stability_window = 3L) {
  group_counts <- as.integer(group_counts)
  clustering_error <- as.numeric(clustering_error)

  if (length(group_counts) != length(clustering_error) || length(group_counts) == 0) {
    stop("group_counts and clustering_error must have the same positive length.")
  }
  if (any(!is.finite(clustering_error)) || any(clustering_error < 0)) {
    stop("clustering_error must contain finite non-negative values.")
  }
  if (is.unsorted(group_counts, strictly = TRUE)) {
    stop("group_counts must be strictly increasing.")
  }
  if (length(stability_threshold) != 1 || !is.finite(stability_threshold) ||
      stability_threshold < 0 || stability_threshold > 1) {
    stop("stability_threshold must be one finite number between 0 and 1.")
  }
  if (length(stability_window) != 1 || !is.finite(stability_window) ||
      stability_window < 1 || stability_window != floor(stability_window)) {
    stop("stability_window must be one positive integer.")
  }

  monotone_error <- cummin(clustering_error)
  relative_improvement <- rep(NA_real_, length(monotone_error))
  if (length(monotone_error) > 1) {
    previous_error <- monotone_error[-length(monotone_error)]
    # Entry K is the proportional error reduction from K - 1 to K groups.
    relative_improvement[-1] <- ifelse(
      previous_error > 0,
      (previous_error - monotone_error[-1]) / previous_error,
      0
    )
  }

  path_length <- length(group_counts)
  effective_window <- min(as.integer(stability_window), max(0L, path_length - 1L))
  max_future_improvement <- rep(NA_real_, path_length)
  stable <- rep(FALSE, path_length)
  error_drop <- monotone_error[1] - monotone_error[path_length]
  error_scale <- max(1, abs(monotone_error[1]))

  if (path_length == 1L ||
      error_drop <= sqrt(.Machine$double.eps) * error_scale) {
    selected_index <- 1L
    stable[1] <- TRUE
    max_future_improvement[1] <- 0
  } else {
    last_candidate <- path_length - effective_window
    for (candidate_index in seq_len(last_candidate)) {
      future_indices <- (candidate_index + 1L):(candidate_index + effective_window)
      max_future_improvement[candidate_index] <- max(
        relative_improvement[future_indices],
        na.rm = TRUE
      )
      stable[candidate_index] <- all(
        relative_improvement[future_indices] <= stability_threshold
      )
    }
    selected_index <- if (any(stable)) which(stable)[1] else path_length
  }

  path <- data.frame(
    groups = group_counts,
    clustering_error = clustering_error,
    monotone_error = monotone_error,
    relative_improvement = relative_improvement,
    max_future_improvement = max_future_improvement,
    stability_threshold = rep(stability_threshold, path_length),
    stability_window = rep(effective_window, path_length),
    stable = stable,
    selected = seq_along(group_counts) == selected_index,
    stringsAsFactors = FALSE
  )

  list(groups = group_counts[selected_index], path = path)
}

# Minimize a scale-free clustering loss plus the cost of the extra parameters
# that the selected groups create in the downstream additive GFE model.
select_penalized_groups <- function(group_counts, clustering_error,
                                    parameter_increment, sample_size,
                                    penalty_multiplier = 1) {
  group_counts <- as.integer(group_counts)
  clustering_error <- as.numeric(clustering_error)

  if (length(group_counts) != length(clustering_error) || length(group_counts) == 0) {
    stop("group_counts and clustering_error must have the same positive length.")
  }
  if (any(!is.finite(clustering_error)) || any(clustering_error < 0)) {
    stop("clustering_error must contain finite non-negative values.")
  }
  if (is.unsorted(group_counts, strictly = TRUE) || any(group_counts < 1L)) {
    stop("group_counts must contain strictly increasing positive integers.")
  }
  if (length(parameter_increment) != 1 || !is.finite(parameter_increment) ||
      parameter_increment <= 0) {
    stop("parameter_increment must be one finite positive number.")
  }
  if (length(sample_size) != 1 || !is.finite(sample_size) || sample_size <= 1) {
    stop("sample_size must be one finite number greater than one.")
  }
  if (length(penalty_multiplier) != 1 || !is.finite(penalty_multiplier) ||
      penalty_multiplier < 0) {
    stop("penalty_multiplier must be one finite non-negative number.")
  }

  monotone_error <- cummin(clustering_error)
  baseline_error <- monotone_error[1]
  path_length <- length(group_counts)
  added_parameters <- parameter_increment * (group_counts - 1L)
  penalty_weight <- penalty_multiplier * log(sample_size)
  complexity_penalty <- penalty_weight * added_parameters / sample_size

  if (path_length == 1L || baseline_error <= .Machine$double.xmin) {
    relative_error <- rep(1, path_length)
    log_relative_error <- rep(0, path_length)
    selected_index <- 1L
  } else {
    error_floor <- max(.Machine$double.xmin, baseline_error * .Machine$double.eps)
    relative_error <- monotone_error / baseline_error
    log_relative_error <- log(pmax(monotone_error, error_floor) / baseline_error)
    selection_score <- log_relative_error + complexity_penalty
    selected_index <- which.min(selection_score)
  }

  selection_score <- log_relative_error + complexity_penalty
  relative_improvement <- rep(NA_real_, path_length)
  if (path_length > 1L) {
    previous_error <- monotone_error[-path_length]
    relative_improvement[-1] <- ifelse(
      previous_error > 0,
      (previous_error - monotone_error[-1]) / previous_error,
      0
    )
  }

  path <- data.frame(
    groups = group_counts,
    clustering_error = clustering_error,
    monotone_error = monotone_error,
    relative_error = relative_error,
    log_relative_error = log_relative_error,
    relative_improvement = relative_improvement,
    added_parameters = added_parameters,
    penalty_weight = rep(penalty_weight, path_length),
    complexity_penalty = complexity_penalty,
    selection_score = selection_score,
    selected = seq_along(group_counts) == selected_index,
    stringsAsFactors = FALSE
  )

  list(groups = group_counts[selected_index], path = path)
}

cluster_general <- function(Y, X_list, N, T, init = 100, type = "long", groups = NULL, cc = 0, gamma = 1, cluster_type = 'kmeans', group_selection = 'variance', max_groups = NULL, stability_threshold = 0.05, stability_window = 3L, penalty_parameter_increment = NULL, penalty_sample_size = NULL, penalty_multiplier = 1) {
  dimtheta <- length(X_list)
  mdim <- dimtheta + 1
  group_selection <- match.arg(group_selection, c('variance', 'elbow', 'stability', 'penalized'))
  cluster_type <- match.arg(cluster_type, c('kmeans', 'kcenter'))
  if (is.null(groups) && group_selection == 'elbow' && cluster_type != 'kcenter') {
    stop("group_selection = 'elbow' requires cluster_type = 'kcenter', because the elbow path is defined from the k-center covering radius.")
  }

  if (type == "long") {
    ## --- Unit-side moments ---
    mom_i <- rowMeans(Y)
    for (X in X_list) mom_i <- cbind(mom_i, rowMeans(X))

    mom_micro <- Y
    for (X in X_list) mom_micro <- cbind(mom_micro, X)

    ## --- Demean & rescale ---
    for (j in 1:mdim) {
      col_mean <- mean(mom_i[, j])
      col_sd <- sd(mom_i[, j])
      if (col_sd > 0) {
        cols <- ((j - 1) * T + 1):(j * T)
        mom_micro[, cols] <- (mom_micro[, cols] - col_mean) / col_sd
        mom_i[, j] <- (mom_i[, j] - col_mean) / col_sd
      }
    }

    ## --- Reweight moments only when requested ---
    if (cc == 1) {
      var_w <- var_tot <- matrix(0, mdim, mdim)
      for (j in seq_len(mdim)) {
        for (k in seq_len(mdim)) {
          cols_j <- ((j - 1) * T + 1):(j * T)
          cols_k <- ((k - 1) * T + 1):(k * T)
          var_w[j, k] <- mean(colMeans((mom_micro[, cols_j] - mom_i[, j]) *
                                         (mom_micro[, cols_k] - mom_i[, k]))) / T
          var_tot[j, k] <- mean(mom_i[, j] * mom_i[, k])
        }
      }
      var_b <- var_tot - var_w
      mat_scali <- diag(diag(var_b)) / diag(diag(var_tot))
      mat_scali <- diag(diag(mat_scali))
      mat_scali <- pmax(mat_scali, 0)

      # Safety check: replace with identity if NaN, Inf, or all zeros
      if (any(is.nan(mat_scali)) || any(is.infinite(mat_scali)) || max(mat_scali) == 0) {
        mat_scali <- diag(mdim)
      }
    } else {
      mat_scali <- diag(mdim)
    }

    # --- Rescale data ---
    moment_weights <- diag(mat_scali)
    mom_i <- sweep(mom_i, 2, moment_weights, `*`)
    mom_micro <- sweep(mom_micro, 2, rep(moment_weights, each = T), `*`)
    dimnames(mom_i) <- NULL
    dimnames(mom_micro) <- NULL

    # --- Variance / noise on rescaled data ---
    expanded_mom_i <- mom_i[, rep(seq_len(mdim), each = T), drop = FALSE]
    variance <- sum((mom_micro - expanded_mom_i)^2) / (N * T^2)
    # sum(sapply(1:T, function(t) {
    #   norm(mom_micro[, ((1:mdim) - 1) * T + t] - mom_i, type = "F")^2
    # })) / (N * T^2)

    data <- mom_i
    dim_size <- N

  } else if (type == "tall") {
    ## --- Time-side moments ---
    mom_t <- colMeans(Y)
    for (X in X_list) mom_t <- cbind(mom_t, colMeans(X))

    mom_micro2 <- t(Y)
    for (X in X_list) mom_micro2 <- cbind(mom_micro2, t(X))

    ## --- Demean & rescale ---
    for (j in 1:mdim) {
      col_mean <- mean(mom_t[, j])
      col_sd <- sd(mom_t[, j])
      if (col_sd > 0) {
        cols <- ((j - 1) * N + 1):(j * N)
        mom_micro2[, cols] <- (mom_micro2[, cols] - col_mean) / col_sd
        mom_t[, j] <- (mom_t[, j] - col_mean) / col_sd
      }
    }

    ## --- Reweight moments only when requested ---
    if (cc == 1) {
      var_w2 <- var_tot2 <- matrix(0, mdim, mdim)
      for (j in seq_len(mdim)) {
        for (k in seq_len(mdim)) {
          cols_j <- ((j - 1) * N + 1):(j * N)
          cols_k <- ((k - 1) * N + 1):(k * N)
          var_w2[j, k] <- mean((mom_micro2[, cols_j] - mom_t[, j]) *
                                 (mom_micro2[, cols_k] - mom_t[, k])) / N
          var_tot2[j, k] <- mean(mom_t[, j] * mom_t[, k])
        }
      }
      var_b2 <- var_tot2 - var_w2
      mat_scali2 <- diag(diag(var_b2)) / diag(diag(var_tot2))
      mat_scali2 <- diag(diag(mat_scali2))
      mat_scali2 <- pmax(mat_scali2, 0)
    } else {
      mat_scali2 <- diag(mdim)
    }

    # Safety check: replace with identity if NaN, Inf, or zeros
    if (any(is.nan(mat_scali2)) || any(is.infinite(mat_scali2)) || max(mat_scali2) == 0) {
      mat_scali2 <- diag(mdim)
    }

    # --- Rescale data ---
    moment_weights <- diag(mat_scali2)
    mom_t <- sweep(mom_t, 2, moment_weights, `*`)
    mom_micro2 <- sweep(mom_micro2, 2, rep(moment_weights, each = N), `*`)
    dimnames(mom_t) <- NULL
    dimnames(mom_micro2) <- NULL


    ## --- Variance / noise on rescaled data ---
    expanded_mom_t <- mom_t[, rep(seq_len(mdim), each = N), drop = FALSE]
    variance <- sum((mom_micro2 - expanded_mom_t)^2) / (T * N^2)
    # sum(sapply(1:N, function(i) {
    #   norm(mom_micro2[, ((1:mdim)-1)*N + i] - mom_t, type = "F")^2
    # })) / (T * N^2)

    data <- mom_t
    dim_size <- T

  } else {
    stop("Invalid type. Use 'long' or 'tall'.")
  }

  ## --- Clustering and group-count selection ---
  data <- as.matrix(data)
  data_sd <- apply(data, 2, sd)
  has_variation <- any(is.finite(data_sd) & data_sd > 0)
  distinct_rows <- nrow(unique(data))
  maximum_available_groups <- max(1L, min(dim_size, distinct_rows))
  theoretical_max_groups <- max(1L, as.integer(floor(2 * dim_size^(1 / 3))))
  automatic_max_groups <- min(maximum_available_groups, theoretical_max_groups)
  requested_max_groups <- NA_integer_
  dmat <- if (cluster_type == 'kcenter') as.matrix(dist(data)) else NULL
  selection_path <- NULL

  fit_groups <- function(group_count) {
    group_count <- max(1L, min(as.integer(group_count), maximum_available_groups))
    if (cluster_type == 'kmeans') {
      stats::kmeans(
        data,
        centers = group_count,
        algorithm = "Lloyd",
        nstart = init,
        iter.max = 100
      )
    } else {
      kcenter(data, centers = group_count, nstart = init, dmat = dmat)
    }
  }

  selection_error <- function(result) {
    if (group_selection == 'elbow') {
      return(as.numeric(result$kcenter_radius))
    }
    as.numeric(result$tot.withinss / dim_size)
  }
  selection_objective <- if (group_selection == 'elbow') {
    'kcenter_radius'
  } else {
    'average_within_cluster_sum_of_squares'
  }

  apply_group_selection <- function(candidate_groups, clustering_error) {
    switch(
      group_selection,
      elbow = select_elbow_groups(candidate_groups, clustering_error),
      stability = select_stability_groups(
        candidate_groups,
        clustering_error,
        stability_threshold = stability_threshold,
        stability_window = stability_window
      ),
      penalized = select_penalized_groups(
        candidate_groups,
        clustering_error,
        parameter_increment = if (is.null(penalty_parameter_increment)) {
          if (type == 'long') T else N
        } else penalty_parameter_increment,
        sample_size = if (is.null(penalty_sample_size)) N * T else penalty_sample_size,
        penalty_multiplier = penalty_multiplier
      )
    )
  }

  label_elbow_path <- function(path) {
    if (is.null(path) || group_selection != 'elbow') {
      return(path)
    }

    if (type == 'long') {
      path$K <- path$groups
      path$R_k <- path$clustering_error
      path$R_k_monotone <- path$monotone_error
      path$x_k <- path$normalized_group_count
      path$y_k <- path$normalized_clustering_error
      path$S_k <- path$elbow_score
    } else {
      path$L <- path$groups
      path$R_l <- path$clustering_error
      path$R_l_monotone <- path$monotone_error
      path$x_l <- path$normalized_group_count
      path$y_l <- path$normalized_clustering_error
      path$S_l <- path$elbow_score
    }

    path
  }

  if (!is.null(groups)) {
    if (length(groups) != 1 || !is.finite(groups) || groups < 1) {
      stop("groups must be one finite positive number.")
    }
    k_result <- fit_groups(floor(groups))
    clusters <- nrow(as.matrix(k_result$centers))
  } else if (!has_variation || maximum_available_groups == 1L) {
    clusters <- 1L
    k_result <- fit_groups(clusters)
    if (group_selection %in% c('elbow', 'stability', 'penalized')) {
      selection_path <- apply_group_selection(
        candidate_groups = clusters,
        clustering_error = selection_error(k_result)
      )$path
      selection_path$clustering_objective <- selection_objective
      selection_path <- label_elbow_path(selection_path)
    }
  } else if (group_selection %in% c('elbow', 'stability', 'penalized')) {
    if (is.null(max_groups)) {
      max_groups <- max(10L, ceiling(sqrt(dim_size)))
    }
    if (length(max_groups) != 1 || !is.finite(max_groups) || max_groups < 2) {
      stop("max_groups must be one finite number greater than or equal to 2.")
    }
    requested_max_groups <- as.integer(floor(max_groups))
    automatic_max_groups <- min(automatic_max_groups, requested_max_groups)
    candidate_groups <- seq_len(automatic_max_groups)
    candidate_results <- lapply(candidate_groups, fit_groups)
    clustering_error <- vapply(
      candidate_results,
      selection_error,
      numeric(1)
    )
    selection <- apply_group_selection(candidate_groups, clustering_error)
    selected_index <- match(selection$groups, candidate_groups)
    k_result <- candidate_results[[selected_index]]
    clusters <- nrow(as.matrix(k_result$centers))
    selection_path <- selection$path
    selection_path$effective_groups <- vapply(
      candidate_results,
      function(result) nrow(as.matrix(result$centers)),
      integer(1)
    )
    selection_path$clustering_objective <- selection_objective
    selection_path <- label_elbow_path(selection_path)
  } else {
    if (!is.null(max_groups)) {
      if (length(max_groups) != 1 || !is.finite(max_groups) || max_groups < 1) {
        stop("max_groups must be one finite positive number.")
      }
      requested_max_groups <- as.integer(floor(max_groups))
      automatic_max_groups <- min(automatic_max_groups, requested_max_groups)
    }
    clusters <- 1L
    repeat {
      k_result <- fit_groups(clusters)
      xx <- k_result$tot.withinss / dim_size
      if (xx <= gamma * variance || clusters >= automatic_max_groups) {
        break
      }
      clusters <- clusters + 1L
    }
    clusters <- nrow(as.matrix(k_result$centers))
  }

  if (!is.null(selection_path)) {
    selection_path$theoretical_max_groups <- theoretical_max_groups
    selection_path$automatic_max_groups <- automatic_max_groups
    selection_path$requested_max_groups <- requested_max_groups
  }

  ## --- Return ---
  list(
    res = k_result,
    clusters = clusters,
    centers = k_result$centers,
    selection_rule = if (is.null(groups)) group_selection else 'fixed',
    selection_path = selection_path,
    selection_objective = if (is.null(selection_path)) NULL else selection_objective,
    theoretical_max_groups = theoretical_max_groups,
    automatic_max_groups = automatic_max_groups,
    requested_max_groups = requested_max_groups,
    noise_variance = variance,
    selected_clustering_error = as.numeric(k_result$tot.withinss / dim_size),
    variance_threshold = gamma * variance,
    variance_constraint_met = as.numeric(k_result$tot.withinss / dim_size) <= gamma * variance
  )
}

custom_round <- function(x) floor(x) + as.integer((x - floor(x)) >= 0.5)

kcenter <- function(data, centers, nstart = 100, dmat = NULL) {
  # Ensure data is a numeric matrix
  data <- as.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  centers <- as.integer(centers)

  if (n == 0) {
    stop("kcenter requires at least one observation.")
  }
  if (centers < 1) {
    stop("kcenter requires centers >= 1.")
  }
  centers <- min(centers, n)

  if (is.null(dmat)) {
    dmat <- as.matrix(dist(data))
  }

  best_result <- NULL
  best_radius <- Inf
  nstart <- max(1, as.integer(nstart))
  start_indices <- if (nstart >= n) {
    seq_len(n)
  } else {
    sample.int(n, nstart)
  }
  data_centered <- sweep(data, 2, colMeans(data), check.margin = FALSE)
  totss <- sum(data_centered^2)

  for(first_center in start_indices){
    # Step 2: Initialize centers (first point randomly)
    kcenters <- integer(centers)
    kcenters[1] <- first_center
    min_dist <- dmat[, first_center]
    min_dist[first_center] <- 0

    # Step 3: Greedy selection of remaining centers
    if(centers > 1){
      for(i in 2:centers){
        candidate_dist <- min_dist
        candidate_dist[kcenters[1:(i-1)]] <- -Inf
        kcenters[i] <- which.max(candidate_dist)
        min_dist <- pmin(min_dist, dmat[, kcenters[i]])
      }
    }

    # Step 4: Assign each point to nearest center
    dist_to_centers <- dmat[, kcenters, drop=FALSE]
    cluster_assignments <- max.col(-dist_to_centers, ties.method = "first")

    # Step 5: Maximum distance to nearest center (k-center radius)
    radius <- max(min_dist)

    # Keep the best start (smallest radius)
    if(radius < best_radius){
      best_radius <- radius

      used_clusters <- sort(unique(cluster_assignments))
      if (length(used_clusters) < centers) {
        cluster_map <- match(cluster_assignments, used_clusters)
        kcenters <- kcenters[used_clusters]
        cluster_assignments <- cluster_map
        centers_eff <- length(used_clusters)
      } else {
        centers_eff <- centers
      }

      # Compute within-cluster sum of squares
      withinss <- sapply(1:centers_eff, function(i){
        pts <- data[cluster_assignments == i, , drop=FALSE]
        if(nrow(pts) == 0) return(0)
        sum(rowSums((pts - matrix(data[kcenters[i], ], nrow=nrow(pts), ncol=d, byrow=TRUE))^2))
      })

      # Cluster sizes
      size <- as.vector(table(factor(cluster_assignments, levels=1:centers_eff)))

      # Between-cluster sum of squares
      overall_center <- colMeans(data)
      centers_mat <- as.matrix(data[kcenters, , drop=FALSE])
      betweenss <- sum(
        rowSums((centers_mat - matrix(overall_center, nrow=centers_eff, ncol=d, byrow=TRUE))^2) * size
      )

      # Save result
      best_result <- list(
        centers = centers_mat,
        cluster = cluster_assignments,
        size = size,
        totss = totss,
        withinss = withinss,
        tot.withinss = sum(withinss),
        betweenss = betweenss,
        iter = 1,
        kcenter_radius = radius
      )
    }
  }

  return(best_result)
}
