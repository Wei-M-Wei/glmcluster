# Data Generation Function
generate <- function(N, T, DGP, rho = 0, kappa = 0, theta = -0.5 , link, nu = 5 ) {
  make_ar1 <- function(N, T, kappa) {
    out <- matrix(0, N, T)
    out[, 1] <- rnorm(N)
    for (t in 2:T) {
      out[, t] <- kappa * out[, t - 1] + sqrt(1 - kappa^2) * rnorm(N)
    }
    out
  }

  # Standard normal truncated to [0, 4]
  rtruncnorm_04 <- function(n) {
    qnorm(runif(n, min = pnorm(0), max = pnorm(10)))
  }

  # Truncated-standard-normal marginals with serial dependence governed by rho
  make_truncnorm_ar1 <- function(T, rho) {
    z <- numeric(T)
    z[1] <- rnorm(1)
    for (t in 2:T) {
      z[t] <- rho * z[t - 1] + sqrt(1 - rho^2) * rnorm(1)
    }
    qnorm(pnorm(0) + pnorm(z) * (pnorm(10) - pnorm(0)))
  }
  E <- matrix(0, N, T)
  U <- matrix(0, N, T)
  E[, 1] <- rnorm(N)
  U[, 1] <- rnorm(N)

  for (t in 2:T) {
    E[, t] <- kappa * E[, t - 1] + sqrt(1 - kappa^2) * rnorm(N)
    U[, t] <- kappa * U[, t - 1] + sqrt(1 - kappa^2) * rnorm(N)
  }

  A <- rtruncnorm_04(N)
  B <- make_truncnorm_ar1(T, rho)

  A_rep <- matrix(rep(A, T), N, T)
  B_rep <- t(matrix(rep(B, N), T, N))


  if (DGP == 1) {
    F <- ( 0.5 * A_rep^(10)  +  0.5 * B_rep^(10) )^(1 / (5) )
    H <- ( 0.5 * A_rep^(10)  +  0.5 * B_rep^(10) )^(1 / (5) )
  } else if (DGP == 2) {
    F <- A_rep^2 + (B_rep * A_rep) + 0.5 * sin(B_rep * A_rep)
    H <- B_rep^2 + (B_rep * A_rep) + 0.5 * sin(B_rep * A_rep)
    # F <- A_rep^2 + (B_rep * A_rep)
    # H <- B_rep^2 + 2*(B_rep * A_rep)
  } else if (DGP == 3) {
    A2 <- matrix(rtruncnorm_04(N * 2), nrow = N, ncol = 2)
    B2 <- matrix(0, nrow = T, ncol = 2)
    for (d in 1:2) {
      B2[, d] <- make_truncnorm_ar1(T, rho)
    }

    A_rep1 <- matrix(rep(A2[, 1], T), N, T)
    A_rep2 <- matrix(rep(A2[, 2], T), N, T)
    B_rep1 <- t(matrix(rep(B2[, 1], N), T, N))
    B_rep2 <- t(matrix(rep(B2[, 2], N), T, N))

    F1 <- ( 0.5 * A_rep1^(10)  +   0.5 * B_rep1^(10) )^(1 / (5) )
    F2 <- (  0.5 * A_rep2^(10) +  0.5 * B_rep2^(10) )^(1 / (5) )
    F <- F1 + F2

    H1 <- (   A_rep1^(10)  +   B_rep1^(10) )^(1 / (5) )
    H2 <- sin(A_rep2 * B_rep2)

    X <- H1 + H2 + U


    linear_predictor <- theta * X + F

    if (link == 'probit'){
      Y <- (linear_predictor + matrix(rnorm(N * T, 0, 1), N, T) > 0)
    }else{
      Y <- (linear_predictor + matrix(rlogis(N * T, 0, 1), N, T) > 0)
    }
    Y <- matrix(as.numeric(Y), nrow = N, ncol = T)
  }else if (DGP == 4) {
    A2 <- matrix(rtruncnorm_04(N * 2), nrow = N, ncol = 2)
    B2 <- matrix(0, nrow = T, ncol = 2)
    for (d in 1:2) {
      B2[, d] <- make_truncnorm_ar1(T, rho)
    }

    A_rep1 <- matrix(rep(A2[, 1], T), N, T)
    A_rep2 <- matrix(rep(A2[, 2], T), N, T)
    B_rep1 <- t(matrix(rep(B2[, 1], N), T, N))
    B_rep2 <- t(matrix(rep(B2[, 2], N), T, N))

    F1 <- (0.5 * A_rep1^(10)  + 0.5 * B_rep1^(10) )^(1 / (5) )
    F2 <- (0.5 * A_rep2^(10)  + 0.5 * B_rep2^(10) )^(1 / (5) )
    F <- F1 + F2

    H1 <- ( A_rep1^(10)  +  B_rep1^(10) )^(1 / (5) )
    H2 <- sin(A_rep2 * B_rep2)

    X1 <- H1 + H2 + U


    H_prime1 <- (0.5 * A_rep1^(10)  +  0.5 * B_rep1^(10) )^(1 / (10) ) # A_rep1^2 + (B_rep1 * A_rep1) + sin(B_rep1 * A_rep1)
    H_prime2 <- (0.5 * A_rep2^(10)  +  0.5 * B_rep2^(10) )^(1 / (10) ) # A_rep2^2 + (B_rep2 * A_rep2) + sin(B_rep2 * A_rep2)
    U_prime <- make_ar1(N, T, kappa)
    X2 <- H_prime1 + H_prime2 + U_prime

    linear_predictor <- theta[1] * X1 + theta[2] * X2 + F

    if (link == 'probit'){
      Y <- (linear_predictor + matrix(rnorm(N * T, 0, 1), N, T) > 0)
    }else{
      Y <- (linear_predictor + matrix(rlogis(N * T, 0, 1), N, T) > 0)
    }
    Y <- matrix(as.numeric(Y), nrow = N, ncol = T)
  }

  if (DGP != 3 && DGP != 4) {
    X <- H + U
    if (link == 'probit'){
      Y <- (X * theta + F + matrix(rnorm(N * T, 0, 1), N, T) > 0)
    }else{
      Y <- (X * theta + F + matrix(rlogis(N * T, 0, 1), N, T) > 0)
    }
    Y <- matrix(as.numeric(Y), nrow = N, ncol = T)
  }

  id <- numeric(N * T)
  for (i in 1:N) {
    for (t in 1:T) {
      id[((t - 1) * N + i):(t * N)] <- i
    }
  }

  time <- numeric(N * T)
  for (i in 1:N) {
    for (t in 1:T) {
      time[((t - 1) * N + 1):(t * N)] <- t
    }
  }

  vY <- c(Y)
  if (DGP == 4) {
    vX1 <- c(X1)
    vX2 <- c(X2)
    data.frame(id, time, vY, vX1, vX2)
  } else {
    vX <- c(X)
    data.frame(id, time, vY, vX)
  }


}
