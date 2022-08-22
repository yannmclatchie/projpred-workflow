# load necessary libraries
library(tidyverse)
library(brms)
library(projpred)
library(parallel)

# load data generating scripts
source("./R/generate.R")

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
  
  # perform projection predictive inference
  cv_loo <- cv_varsel(
    fit, 
    method = "forward", 
    cv_method = "LOO",
    validate_search = TRUE, 
    nterms_max = 25
  )
  cv_kfold <- cv_varsel(
    fit, 
    method = "forward", 
    cv_method = "kfold",
    validate_search = TRUE, 
    nterms_max = 25
  )

  # extract projection summaries
  loo_df <- summary(cv_loo, stats = c("elpd", "rmse"))$selection
  loo_df["method"] = "forward"
  loo_df["loo_wts"] <- exp(loo_df["elpd"]) / sum(exp(loo_df["elpd"]))
  kfold_df <- summary(cv_kfold, stats = c("elpd", "rmse"))$selection

  # compute LOO weights across submodels
  summary["loo_wts"] <- exp(summary["elpd"]) / sum(exp(summary["elpd"]))

  # write the table
  csv_name <- paste0("./data/elbows/",n,"n,",rho,"rho.csv")
  write.table(summary, file = csv_name, sep = ",", row.names = FALSE)
}

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
