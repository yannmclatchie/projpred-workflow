# load necessary libraries
library(dplyr)
library(brms)
library(projpred)
library(parallel)

# load data generating scripts
source("./R/generate.R")

# simulate data
rho = 0.5
n = 200
data <- .generate(rho, n)
indep_data <- .generate(rho, n)
d_test <- list(
  "data" = indep_data,
  "offset" = rep(0, n),
  "weights" = rep(1, n),
  "y" = indep_data[, names(indep_data) %in% c("y")]
)

# fit the BRMS model
fit <- brm(
  formula = "y ~ .", 
  family = "gaussian", 
  prior = set_prior("horseshoe(11)"), 
  data = data
)

# perform projection predictive inference
cv_vs <- cv_varsel(
  fit, 
  d_test = d_test,
  method = "forward", 
  nterms_max = 25,
  cv_method = "kfold",
  K = 10,
  validate_search = TRUE
)
plot(cv_vs, stat="mlpd")

# define experiment
experiment <- function(rho, n){
  
  # simulate data
  data <- .generate(rho, n)
  indep_data <- .generate(rho, n)
  d_test <- list(
    "data" = indep_data,
    "offset" = rep(0, n),
    "weights" = rep(1, n),
    "y" = indep_data[, names(indep_data) %in% c("y")]
  )
  
  # fit the BRMS model
  fit <- brm(
    formula = "y ~ .", 
    family = "gaussian", 
    prior = set_prior("horseshoe(11)"), 
    data = data
  )
  
  # perform projection predictive inference
  vs <- varsel(
    fit, 
    d_test = d_test,
    method = "forward", 
    nterms_max = 25,
    cv_method = "kfold",
    k = 10
  )
  
  # compute submodel summaries
  df <- summary(vs, stats = c("elpd", "rmse"))$selection
  df["rho"] <- rho
  df["n"] <- n
  df["loo_wts"] <- exp(df["elpd"]) / sum(exp(df["elpd"]))
  
  # write the table
  csv_name <- paste0("./data/elbows/",n,"n,",rho,"rho.csv")
  ff <- file(csv_name, open="w")
  write.table(df, file = ff, sep = ",", row.names = FALSE)
  close(ff)
}

# vary both the number of data observations and the correlation between parameters
rhos <- c(0, 0.5, 0.9)
ns <- c(50, 100, 200, 500)
options <- expand.grid(rho=rhos, n=ns)
