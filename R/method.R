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
  vs_forward <- varsel(
    fit, 
    d_test = d_test,
    method = "forward", 
    nterms_max = 25
  )
  vs_l1 <- varsel(
    fit, 
    d_test = d_test,
    method = "L1", 
    nterms_max = 25
  )
  
  # compute reference model criteria on independent test data
  df_forward <- summary(vs_forward, stats = c("elpd", "rmse"))$selection
  df_forward["rho"] = rho
  df_forward["n"] = n
  df_forward["method"] = "forward"
  df_l1 <- summary(vs_l1, stats = c("elpd", "rmse"))$selection
  df_l1["rho"] = rho
  df_l1["n"] = n
  df_l1["method"] = "L1"

  # initialise summary dataframe
  summary <- rbind(df_forward, df_l1)
  
  # write the table
  csv_name <- paste0("./data/method/regression/",n,"n,",rho,"rho.csv")
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
