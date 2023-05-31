# load necessary libraries
library(Matrix)
library(MASS)

# convert a correlation matrix to a covariance matrix
.cor2cov <- function(C, var = NULL) {
  if (is.null(var)) stop(
    "cor2cov: cannot calculate covariance matrix without variances"
  )
  if (ncol(C) != nrow(C)) stop(
    "cor2cov: 'C' is not a square numeric matrix!"
  )
  if (length(var) != ncol(C)) stop(
    "cor2cov: length of 'var' and dimension of 'C' are not equal!"
  )
  if (any(!is.finite(var))) warning(
    "cor2cov: 'var' had 0 or NA entries; result is doubtful!"
  )
  d <- sqrt(var)
  V <- outer(d, d) * C
  return(V)
}

# define the data generating function
.generate <- function(rho, n) {

  # number of parameters in the model
  p <- 100

  # define xi given rho
  if (rho == 0){
    xi <- 0.59
  } else if (rho == 0.5){
    xi <- 0.34
  } else if (rho == 0.9){
    xi <- 0.28
  } else {
    warning("Invalid rho provided.")
  }

  # build correlation matrix
  block_size <- 5
  num_matrices <- p / block_size
  listOfMatrices <- vector("list", num_matrices)
  for (i in 1:num_matrices) {
    listOfMatrices[[i]] <- matrix(
      rep(rho, block_size * block_size), nrow=block_size
    )
  }
  R <- matrix(bdiag(listOfMatrices), nrow=p)
  diag(R) <- 1

  # build associated covariate weights
  w1 <- rep(xi, block_size)
  w2 <- rep(xi * 0.5, block_size)
  w3 <- rep(xi * 0.25, block_size)
  w4 <- rep(0, (num_matrices - 3) * block_size)
  w <- c(w1, w2, w3, w4)

  # define zero mean vector
  sigma <- 1
  mu <- rep(0, p)
  sds <- rep(sigma, p)
  S <- .cor2cov(R, sds)

  # sample data
  X <- mvrnorm(n=n, mu=mu, Sigma=S)
  y <- X %*% w + sigma * rnorm(n=n)
  data <- data.frame(y = y, X = X)

  return(data)
}

.binom_samples <- function (n, p) {
  return( rbinom(1, n, p) )
}

.invlogit <- function(x) {
  # Taken from LaplacesDemon
  InvLogit <- 1 / {1 + exp(-x)}
  return(InvLogit)
}

.generate_binomial <- function (rho, N) {
  # define number of covariates
  p <- 100

  # define xi given rho
  if (rho == 0){
    xi <- 0.59
  } else if (rho == 0.5){
    xi <- 0.34
  } else if (rho == 0.9){
    xi <- 0.28
  } else {
    warning("Invalid rho provided.")
  }

  # build correlation matrix
  block_size <- 5
  num_matrices <- p / block_size
  listOfMatrices <- vector("list", num_matrices)
  for (i in 1:num_matrices) {
    listOfMatrices[[i]] <- matrix(
      rep(rho, block_size * block_size), nrow = block_size
    )
  }
  R <- matrix(Matrix::bdiag(listOfMatrices), nrow = p)
  diag(R) <- 1

  # build associated covariate weights
  w1 <- rep(xi, block_size)
  w2 <- rep(xi * 0.5, block_size)
  w3 <- rep(xi * 0.25, block_size)
  w4 <- rep(0, (num_matrices - 3) * block_size)
  w <- c(w1, w2, w3, w4)

  # define zero mean vector
  sigma <- 1
  mu <- rep(0, p)
  sds <- rep(sigma, p)
  S <- .cor2cov(R, sds)

  # sample covariates and compute latent predictor
  X <- MASS::mvrnorm(n = N, mu = mu, Sigma = S)
  eta <- X %*% w + sigma * rnorm(n = N)

  # apply inverse link function to build probabilities
  ps <- .invlogit(eta)

  # sample from variate
  y <- mapply(.binom_samples, 1, ps)

  data <- data.frame(y = y, X = X)
}

#' Generate data with fixed signal-to-noise ratio (R^2)
#'
#' @param n number of data points
#' @param p number of predictors
#' @param sigma residual variance
#' @param R2 R^2 of the data generating process
#' @param rho correlation between predictors
#' @param beta for convenience, pre-computed beta values can also be provided.
#'     Then, only design-matrix and noise are re-drawn.
#' @param fixed_beta Whether to fix betas into a specific value, so that the DGP
#'     has a R2 of the desired value. If true, the variance of beta vector is
#'     solved and values of beta are drawn randomly.
#' @return list with two elements, `df` and `beta`. The former is data frame
#'     with simulated predictors and response and latter the coefficient values.
generate_fixed_R2_data = function(
    n,
    p,
    sigma,
    R2,
    rho = NULL,
    beta = NULL,
    fixed_beta = FALSE) {
  stopifnot(R2 >= 0, R2 < 1)
  X = generate_X(n, p, rho)
  if(R2 == 0) {
    beta = rep(0, p)
  } else {
    if(fixed_beta) { beta = .fixed_R2_fixed_beta(n, p, sigma, R2, rho) }
    if(!fixed_beta) { beta = .fixed_R2_random_beta(n, p, sigma, R2, beta) }
  }
  eps = rnorm(n, 0, sigma)
  y = X %*% beta + eps
  return(list(df = data.frame(X, y), beta = beta))
}


.fixed_R2_random_beta = function(n, p, sigma, R2, beta) {
  if(is.null(beta)) {
    var_beta = sigma^2/p / (1/R2 - 1)
    beta = rnorm(p, 0, sqrt(var_beta))
  }
  beta
}


.fixed_R2_fixed_beta = function(n, p, sigma, R2, rho) {
  beta = sqrt(sigma^2 * R2 / ((1 - R2)*p*(1 + p*rho - rho)))
  rep(beta, p)
}


generate_X = function(n, p, rho = NULL) {
  mu = rep(0, p)
  Sigma = diag(p)
  if(!is.null(rho)) {
    Sigma[Sigma != 1] = rho
  }

  MASS::mvrnorm(n, mu, Sigma)
}
