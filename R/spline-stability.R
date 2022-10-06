# check if scam if installed
if (!require("scam")) install.packages("scam", repos = "http://cran.us.r-project.org"); library(scam)

# set seed for reproducibility
set.seed(210426)

# load necessary libraries
library(tidyverse)
library(brms)
library(projpred)
library(parallel)
library(gamm4)

# load data generating scripts
source("./R/generate.R")

# set number of repetitions 
num_runs <- 10

# define experiment
experiment <- function(rho, n){
  
  df <- data.frame(matrix(ncol = 18, nrow = 0))
  
  for (run in 1:num_runs) {
    # simulate data
    data <- .generate(rho, n)
    
    # fit the BRMS model with an R2D2 prior
    fit <- brm(
      formula = "y ~ .", 
      family = "gaussian", 
      prior =set_prior(R2D2(mean_R2 = 0.5, prec_R2 = 50, cons_D2 = 1, autoscale = TRUE)), 
      data = data
    )
    
    # perform projection predictive inference
    vs <- varsel(
      fit, 
      method = "forward", 
      nterms_max = 100
    )
    
    # compute submodel summaries
    run_df <- summary(vs, stats = c("elpd"))$selection
    
    # extract reference model ELPD
    ref.loo.elpd <- run_df %>%
      dplyr::arrange(as.numeric(size)) %>%
      dplyr::slice(n = dplyr::n()) %>%
      dplyr::select(elpd)
    run_df["ref.loo.elpd"] <- ref.loo.elpd
    
    # fit different splines to the projpred results
    spline <- gam(elpd/se ~ s(size, k = 25, m = 2), data = run_df)
    mono.spline <- scam(elpd/se ~ s(size, k = 25, bs = "mpi", m = 2), data = run_df)
    spline.diff <- gam(diff/diff.se ~ s(size, k = 25, m = 2), data = run_df)
    mono.spline.diff <- scam(diff/diff.se ~ s(size, k = 25, bs = "mpi", m = 2), data = run_df)
    
    # update dataframe
    run_df <- run_df %>%
      mutate(
        spline.elpd = spline$fit*run_df[,'se'],
        spline.elpd.se = sqrt(spline$sig2)*run_df[,'se'],
        mono.spline.elpd = mono.spline$fit*run_df[,'se'],
        mono.spline.elpd.se = sqrt(mono.spline$sig2)*run_df[,'se'],
        spline.elpd.diff = spline.diff$fit*run_df[,'diff.se'],
        spline.elpd.diff.se = sqrt(spline.diff$sig2)*run_df[,'diff.se'],
        mono.spline.elpd.diff = mono.spline.diff$fit*run_df[,'diff.se'],
        mono.spline.elpd.diff.se = sqrt(mono.spline.diff$sig2)*run_df[,'diff.se']
      )
    run_df["n"] <- n
    run_df["rho"] <- rho
    run_df["run"] <- run
    
    # update results
    df <- rbind(df, run_df)
  }
  return (df)
}

# vary both the number of data observations and the correlation between parameters
rhos <- c(0, 0.5, 0.9)
ns <- c(100, 200, 500)
options <- expand.grid(rho=rhos, n=ns)

# perform experiment
res <- parallel::mcMap(
  f = experiment, 
  rho = options[, "rho"], 
  n = options[, "n"], 
  mc.cores = Sys.getenv('SLURM_CPUS_PER_TASK')
)

# concatenate experiment results
df <- do.call("rbind", res)

# write the table
setwd("/scratch/work/mclatcy1/projpred-workflow")
csv_name <- paste0("./data/smoothed/spline_stability.csv")
ff <- file(csv_name, open="w")
write.table(df, file = ff, sep = ",", row.names = FALSE)
close(ff)
