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
  
  message("Fitting 10 ...")
  # using recommended parameters
  ref_model_10 <- cv_varsel(
    fit, 
    d_test = d_test,
    nterms_max = p,
    ndraws_pred = 10,
    nclusters_pred = NULL,
    method = "forward"
  )
  summary <- summary(ref_model_10, stats = c("elpd", "rmse"))$selection
  summary["ndraws_pred"] = 10
  summary["rho"] = rho
  summary["n"] = n
  
  message("Fitting 200 ...")
  # using fewer pred draws
  ref_model_200 <- cv_varsel(
    fit, 
    d_test = d_test,
    nterms_max = p,
    ndraws_pred = 200,
    nclusters_pred = NULL,
    method = "forward"
  )
  summary200 <- summary(ref_model_200, stats = c("elpd", "rmse"))$selection
  summary200["ndraws_pred"] = 200
  summary200["rho"] = rho
  summary200["n"] = n
  summary <- rbind(summary, summary200)
  
  message("Fitting 50 ...")
  ref_model_50 <- cv_varsel(
    fit, 
    d_test = d_test,
    nterms_max = p,
    ndraws_pred = 50,
    nclusters_pred = NULL,
    method = "forward"
  )
  summary50 <- summary(ref_model_50, stats = c("elpd", "rmse"))$selection
  summary50["ndraws_pred"] = 50
  summary50["rho"] = rho
  summary50["n"] = n
  summary <- rbind(summary, summary50)
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
csv_name <- paste0("./data/pred_draws/regression.csv")
ff <- file(csv_name, open="w")
write.table(df, file = ff, sep = ",", row.names = FALSE)
close(ff)
