generate_X = function(n, p, rho = NULL) {
  mu = rep(0, p)
  Sigma = diag(p)
  if(!is.null(rho)) {
    Sigma[Sigma != 1] = rho
  }
  
  MASS::mvrnorm(n, mu, Sigma)
}


#' Generate data with fixed signal-to-noise ratio (R^2)
#'
#' @param n number of data points
#' @param p number of predictors
#' @param sigma residual variance (if applicable)
#' @param SNR signal-to-noise ratio, or R^2 of the generated data
#' @param rho correlation between predictors
#' @param beta optionally, fixed values of coefficients. If not provided, new 
#'     values are sampled
#' @param family not used atm
#'
#' @return
fixed_SNR = function(
    n, 
    p, 
    sigma, 
    SNR, 
    rho = NULL, 
    beta = NULL, 
    family = 'gaussian') {
  stopifnot(family %in% c('gaussian'),
            SNR >= 0)
  if(family == 'gaussian') return(fixed_SNR_gaussian(n, p, sigma, SNR, rho, beta))
  if(family == 'bernoulli') return(fixed_SNR_bernoulli(n, p, sigma, SNR, rho, beta))
}

fixed_SNR_gaussian = function(n, p, sigma, SNR, rho, beta) {
  X = generate_X(n, p, rho)
  if (is.null(beta)) {
    if(SNR > 0) {
      # Assumes var_x = 1
      var_beta = sigma^2/p / (1/SNR - 1)   
      beta = rnorm(p, 0, sqrt(var_beta))
    } else {
      # If SNR is zero, all coefs are set to zero aswell
      beta = rep(0, p)
    } 
  }
  eps = rnorm(n, 0, sigma)
  y = X %*% beta + eps
  return(list(df = data.frame(X, y), beta = beta))
}


fixed_SNR_bernoulli = function(n, p, sigma, SNR, beta) {
  stop("not implemented yet")
  # X = generate_X(n, p)
  # if (is.null(beta)) {
  #   if(SNR > 0) {
  #     # Assumes var_x = 1 and rho = 0
  #     var_beta = sigma^2/p / (1/SNR - 1)   
  #     beta = rnorm(p, 0, sqrt(var_beta))
  #   } else {
  #     # If SNR is zero, all coefs are set to zero aswell
  #     beta = rep(0, p)
  #   } 
  # }
  # y = rbernoulli(n, 1 / (1 + exp(-X %*% beta)))
  # return(list(df = data.frame(X, y), beta = beta))
}


