# load necessary libraries
library(brms)
library(projpred)

# load data generating scripts
source("./generate.R")

n <- 100
rho <- 0.5
# simulate data
data <- .generate(rho, n)
indep_data <- .generate(rho, n)
d_test <- list(
  # "data" = indep_data[, !names(indep_data) %in% c("y")],
  "data" = indep_data,
  "offset" = rep(0, n),
  "weights" = rep(1, n),
  "y" = indep_data[, names(indep_data) %in% c("y")]
)
  
# define number of parameters in reference model
p <- 25
  
# fit the BRMS model
fit <- brm(
  formula = "y ~ .", 
  family = "gaussian", 
  prior = set_prior("horseshoe(11)"), 
  data = data
)

vs <- varsel(
  fit, 
  d_test = d_test,
  nterms_max = p
)
vs
