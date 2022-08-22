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