# load necessary libraries
library(tidyverse)
library(brms)
library(projpred)
library(parallel)

# load data generating scripts
source("./R/generate.R")

# define experiment
experiment <- function(rho, n){
  
  # fit the BRMS model with an R2D2 prior
  fit <- brm(
    formula = "y ~ .", 
    family = "gaussian", 
    prior = set_prior(R2D2(mean_R2 = 0.25, prec_R2 = 20, cons_D2 = 0.25)), 
    data = data
  )
  
  # perform projection predictive inference
  vs <- varsel(
    fit, 
    method = "forward", 
    nterms_max = 10
  )
  
  summary <- data.frame(matrix(nrow = 0, ncol = 5))
  
  for (size in 0:length(vs$solution_terms)) {
    # define formula
    sel_terms <- head(vs$search_path$solution_terms, size)
    if (size > 0){
      new_formula <- paste0("y ~ ", paste0(sel_terms, collapse = " + "))
    }
    else if (size == 0){
      new_formula <- "y ~ 1"
    }
    else {
      stop("Incorrect size")
    }
    
    # update the brms model
    new_fit <- update(fit, formula. = new_formula)
    
    # compute the log-likelihood on training data
    train_llk <- log_lik(newfit)
    train_loo <- loo::elpd(train_llk)
    train_mlpd <- train_loo$estimates["elpd", "Estimate"] / nrow(fit$data)
    
    # compute the log-likelihood on holdout data
    test_llk <- log_lik(newfit, newdata = newdata)
    test_loo <- loo::elpd(test_llk)
    test_mlpd <- test_loo$estimates["elpd", "Estimate"] / nrow(newfit$data)
    
    # define submodel summary tables
    submdl_train_summary <- data.frame(
      size = size, 
      n = n, 
      rho = rho,
      mlpd = train_mlpd,
      data = "train"
    )
    submdl_test_summary <- data.frame(
      size = size, 
      n = n, 
      rho = rho,
      mlpd = test_mlpd,
      data = "test"
    )
    submdl_summary <- rbind(submdl_train_summary, submdl_test_summary)
    summary <- rbind(submdl_test_summary, summary)
  }
  
  return(summary)
}

# vary both the number of data observations and the correlation between parameters
rhos <- c(0, 0.5, 0.9)
ns <- c(50, 100, 200, 500)
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
csv_name <- paste0("./data/bma/bma_ref.csv")
ff <- file(csv_name, open="w")
write.table(df, file = ff, sep = ",", row.names = FALSE)
close(ff)