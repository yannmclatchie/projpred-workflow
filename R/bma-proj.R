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
  
  # fit the BRMS model with an R2D2 prior
  fit <- brm(
    formula = "y ~ .", 
    family = "gaussian", 
    prior = set_prior(R2D2(mean_R2 = 0.25, prec_R2 = 20, cons_D2 = 0.25)), 
    data = data
  )
  
  # perform projection predictive inference
  vs_insample <- varsel(
    fit, 
    method = "forward", 
    nterms_max = 25
  )
  
  # perform projection predictive inference
  vs_test <- varsel(
    fit, 
    d_test = d_test,
    method = "forward", 
    nterms_max = 25
  )
  
  # compute submodel summaries
  df <- summary(vs_insample, stats = c("mlpd"))$selection
  df["rho"] <- rho
  df["n"] <- n
  df["data"] <- "train"
  
  df_test <- summary(vs_test, stats = c("mlpd"))$selection
  df_test["rho"] <- rho
  df_test["n"] <- n
  df_test["data"] <- "test"
  df <- rbind(df, df_test)
  
  # write the table
  csv_name <- paste0("./data/bma-proj/",n,"n,",rho,"rho.csv")
  ff <- file(csv_name, open="w")
  write.table(df, file = ff, sep = ",", row.names = FALSE)
  close(ff)
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