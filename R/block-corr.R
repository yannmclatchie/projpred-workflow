library(tidyverse)
library(Matrix)
library(MASS)
library(brms)
library(projpred)
library(parallel)

# convert a correlation matrix to a covariance matrix
.cor2cov <- function(C, var = NULL)
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
.generate <- function(rho, n){
  
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
  
  return(data)
}

# define experiment
experiment <- function(rho, n){
  
  # simulate data
  data <- .generate(rho, n)
  
  # fit the BRMS model
  fit <- brm(
    formula = "y ~ .", 
    family = "gaussian", 
    prior = set_prior("horseshoe(11)"), 
    data = data
  )
  
  # simulate new independent test data
  indep_data <- .generate(rho, n)
  
  # build reference model using new independent test data
  ref_model <- init_refmodel(
    fit,
    data = indep_data,
    formula = formula("y ~ ."), 
    family = gaussian(),
    extract_model_data = function(object, newdata = NULL, wrhs = NULL,
                                  orhs = NULL, extract_y = TRUE) {
      if (!extract_y) {
        resp_form <- NULL
      } else {
        resp_form <- ~ y
      }
      
      if (is.null(newdata)) {
        newdata <- indep_data
      }
      
      args <- projpred:::nlist(object, newdata, wrhs, orhs, resp_form)
      return(projpred::do_call(projpred:::.extract_model_data, args))
    },
    cvfun = function(folds) {
      kfold(
        fit, 
        K = max(folds), 
        save_fits = TRUE, 
        folds = folds
      )$fits[, "fit"]
    },
    dis = as.matrix(fit)[, "sigma"]
  )
  
  # using recommended parameters
  cv <- cv_varsel(ref_model, validate_search = TRUE, nterms_max = 25)
  
  # build table of interest
  summary <- summary(cv, stats = c("elpd", "mlpd", "mse", "rmse"))$selection
  summary["n"] <- n
  summary["rho"] <- rho
  summary["kl"] <- cv$kl
  summary["loo_wts"] <- exp(summary["elpd.loo"]) / sum(exp(summary["elpd.loo"]))
  
  # write the table
  csv_name <- paste0("./data/elbows/",n,"n,",rho,"rho.csv")
  write.table(summary, file = csv_name, sep = ",", row.names = FALSE)
}


experiment(0, 50)

# vary both the number of data observations and the correlation between parameters
rhos <- c(0, 0.5, 0.9)
ns <- c(50, 100, 200, 500)
options <- expand.grid(rho=rhos, n=ns)

# perform experiment
parallel::mcmapply(
  FUN = experiment, 
  rho = options[,1], 
  n = options[,2], 
  mc.cores = Sys.getenv('SLURM_CPUS_PER_TASK')
)