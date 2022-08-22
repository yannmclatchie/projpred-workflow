# load necessary libraries
library(tidyverse)
library(brms)
library(projpred)
library(parallel)

# load data generating scripts
source("./R/generate.R")

# define regression experiment
regression_experiment <- function(rho, n){
  
  # simulate data
  data <- .generate(rho, n)
  indep_data <- .generate(rho, n)
  d_test <- list(
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
  
  message("Fitting 0.0 ...")
  # using recommended parameters
  ref_model_0 <- cv_varsel(
    fit, 
    d_test = d_test,
    nterms_max = p,
    method = "L1",
    regul = 1e-04

  )
  summary <- summary(ref_model_0, stats = c("elpd", "rmse"))$selection
  summary["regul"] <- 1e-04
  summary["rho"] = rho
  summary["n"] = n
  
  message("Fitting 0.5 ...")
  ref_model_0.5 <- cv_varsel(
    fit, 
    d_test = d_test,
    nterms_max = p,
    method = "L1",
    regul = 0.5
  )
  summary0.5 <- summary(ref_model_0.5, stats = c("elpd", "rmse"))$selection
  summary0.5["regul"] = 0.5
  summary0.5["rho"] = rho
  summary0.5["n"] = n
  summary <- rbind(summary, summary0.5)
  
  message("Fitting 1 ...")
  ref_model_1 <- cv_varsel(
    fit, 
    d_test = d_test,
    nterms_max = p,
    method = "L1",
    regul = 1 - 1e-04

  )
  summary1 <- summary(ref_model_1, stats = c("elpd", "rmse"))$selection
  summary1["regul"] = 1 - 1e-04
  summary1["rho"] = rho
  summary1["n"] = n
  summary <- rbind(summary, summary1)
  message("Done!")
  
  return (summary)
}

# vary both the number of data observations and the correlation between parameters
rhos <- c(0, 0.5, 0.9)
ns <- c(50, 100, 200, 500)
options <- expand.grid(rho=rhos, n=ns)

# perform experiment
res <- parallel::mcMap(
  f = regression_experiment, 
  rho = options[,1],
  n = options[,2], 
  mc.cores = Sys.getenv('SLURM_CPUS_PER_TASK')
)

# concatenate experiment results
df <- do.call("rbind", res)

# write the table
csv_name <- paste0("./data/regul/regression.csv")
ff <- file(csv_name, open="w")
write.table(df, file = ff, sep = ",", row.names = FALSE)
close(ff)
