# load necessary libraries
library(dplyr)
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
  forward_cv_vs <- cv_varsel(
    fit, 
    d_test = d_test,
    method = "forward", 
    nterms_max = 100
  )
  l1_cv_vs <- cv_varsel(
    fit, 
    d_test = d_test,
    method = "l1", 
    nterms_max = 100
  )
  forward_vs <- varsel(
    fit, 
    d_test = d_test,
    method = "forward", 
    nterms_max = 100
  )
  l1_vs <- varsel(
    fit, 
    d_test = d_test,
    method = "l1", 
    nterms_max = 100
  )
  
  # compute submodel summaries
  forward_cv_summary <- summary(forward_cv_vs, stats = c("elpd", "rmse"))$selection
  l1_cv_summary <- summary(l1_cv_vs, stats = c("elpd", "rmse"))$selection
  forward_summary <- summary(forward_vs, stats = c("elpd", "rmse"))$selection
  l1_summary <- summary(l1_vs, stats = c("elpd", "rmse"))$selection
  
  # define options
  forward_cv_summary["method"] <- "forward"
  forward_summary["method"] <- "forward"
  l1_cv_summary["method"] <- "L1"
  l1_summary["method"] <- "L1"
  forward_cv_summary["validated"] <- "TRUE"
  forward_summary["validated"] <- "FALSE"
  l1_cv_summary["validated"] <- "TRUE"
  l1_summary["validated"] <- "FALSE"
  
  # build and return dataframe
  df <- rbind(forward_cv_summary, forward_summary, l1_cv_summary, l1_summary)
  df["rho"] <- rho
  df["n"] <- n
  return (df)
}

# vary both the number of data observations and the correlation between parameters
rhos <- c(0, 0.5, 0.9)
ns <- c(50, 100, 200, 500)
options <- expand.grid(rho=rhos, n=ns)

# perform experiment
res <- parallel::mcMap(
  f = experiment, 
  rho = options[,1],
  n = options[,2],
  mc.cores = Sys.getenv('SLURM_CPUS_PER_TASK')
)

# concatenate experiment results
df <- do.call("rbind", res)

# save the concatenated results
setwd("/scratch/work/mclatcy1/projpred-workflow")
csv_name <- paste0("./data/stopping.csv")
ff <- file(csv_name, open="w")
write.table(df, file = ff, sep = ",", row.names = FALSE)
close(ff)

