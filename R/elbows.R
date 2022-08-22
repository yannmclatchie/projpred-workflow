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
  vs <- varsel(
    fit, 
    d_test = d_test,
    method = "forward", 
    nterms_max = 25
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

# perform experiment
parallel::mcmapply(
  FUN = experiment, 
  rho = options[,1], 
  n = options[,2], 
  mc.cores = Sys.getenv('SLURM_CPUS_PER_TASK')
)

# define experiment
# experiment_old <- function(rho, n){
  
#   # simulate data
#   data <- .generate(rho, n)
  
#   # fit the BRMS model
#   fit <- brm(
#     formula = "y ~ .", 
#     family = "gaussian", 
#     prior = set_prior("horseshoe(11)"), 
#     data = data
#   )
  
#   # perform projection predictive inference
#   cv <- cv_varsel(
#     fit, 
#     method = "forward", 
#     validate_search = TRUE, 
#     nterms_max = 25
#   )

#   # compute reference model criteria on independent test data
#   indep_data <- .generate(rho, n)
#   new_llk <- log_lik(fit, newdata = newdata)
#   ref.preds <- predict(fit, newdata = newdata)[, "Estimate"]
#   ref.elpd <- loo::elpd(new_llk)$estimates["elpd", "Estimate"]
#   ref.rmse <- sqrt(mean((indep_data$y - ref.preds)^2))

#   # initialise summary dataframe
#   summary <- data.frame(matrix(nrow = 0, ncol = 9))

#   # iterate over solution path
#   for (size in 0:length(cv$solution_terms)) {

#     # project onto subset
#     sel_terms <- head(cv$search_path$solution_terms, size)
#     prj <- project(
#       fit,
#       solution_terms = sel_terms
#     )
    
#     # compute predictions on new independent test data
#     prj_preds <- proj_linpred(prj, newdata = indep_data, integrated = TRUE)
    
#     # compute submodel criteria
#     rmse <- sqrt(mean((indep_data$y - prj_preds$pred)^2))
#     rmse.diff <- ref.rmse - rmse
#     kl <- cv$kl[size + 1]
#     prj.elpd <- loo::elpd(prj_preds$lpd)
#     elpd <- prj.elpd$estimates["elpd", "Estimate"]
#     elpd.se <- prj.elpd$estimates["elpd", "SE"]
#     elpd.diff <- ref.elpd - elpd
    
#     # append to dataframe
#     submdl_summary <- data.frame(
#       size = size, 
#       n = n, 
#       rho = rho,
#       kl = kl,
#       elpd = elpd, 
#       elpd.se = elpd.se, 
#       elpd.diff = elpd.diff,
#       rmse = rmse,
#       rmse.diff = rmse.diff
#     )
#     summary <- rbind(summary, submdl_summary)
#   }

#   # compute LOO weights across submodels
#   summary["loo_wts"] <- exp(summary["elpd"]) / sum(exp(summary["elpd"]))

#   # write the table
#   csv_name <- paste0("./data/elbows/",n,"n,",rho,"rho.csv")
#   write.table(summary, file = csv_name, sep = ",", row.names = FALSE)
# }
