library(tidyverse)
library(Matrix)
library(MASS)
library(parallel)

# convert a correlation matrix to a covariance matrix
cor2cov <- function(C, var = NULL)
{
  if (is.null(var)) stop("cor2cov: cannot calculate covariance matrix without variances")
  if (ncol(C) != nrow(C)) stop("cor2cov: 'C' is not a square numeric matrix!")
  if (length(var) != ncol(C)) stop("cor2cov: length of 'var' and dimension of 'C' are not equal!")
  if (any(!is.finite(var))) warning("cor2cov: 'var' had 0 or NA entries; result is doubtful!")
  d <- sqrt(var)
  V <- outer(d, d) * C
  return(V) 
}

# define the data generating function
generate <- function(rho, n){
  
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
  print(paste("rho = ",rho,". n = ",n))

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
  S <- cor2cov(R, sds)

  # sample datap
  X <- mvrnorm(n=n, mu=mu, Sigma=S)
  y <- X %*% w + sigma * rnorm(n=n)
  data <- data.frame(y = y, X = X)

  # write the data to csv
  csv_name <- paste0("./data/generated/",n,"n,",rho,"rho.csv")
  write.table(data, file = csv_name, sep = ",", row.names = FALSE)
}

# vary both the number of data observations and the correlation between parameters
rhos <- c(0, 0.5, 0.9)
ns <- c(50, 100, 200, 500)
options <- expand.grid(rho=rhos, n=ns)

# perform experiment
parallel::mcmapply(
  FUN = generate, 
  rho = options[,1], 
  n = options[,2], 
  mc.cores = Sys.getenv('SLURM_CPUS_PER_TASK')
)