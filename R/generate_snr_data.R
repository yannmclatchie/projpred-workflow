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