# load necessary libraries
library(tidyverse)
library(brms)
library(projpred)
library(parallel)

# load data generating scripts
source("./R/generate.R")

experiment <- function (rho, N) {
  # simulate data
  data <- .generate_binomial(rho = rho, N = N)
  indep_data <- .generate_binomial(rho = rho, N = N)
  d_test <- list(
    "data" = indep_data,
    "offset" = rep(0, N),
    "weights" = rep(1, N),
    "y" = indep_data[, names(indep_data) %in% c("y")]
  )
  
  # fit binomial family GLM
  fit <- brms::brm(
    formula = "y ~ .", 
    family = bernoulli(link = "logit"),
    prior = brms::set_prior("horseshoe(11)"), 
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
  df <- summary(vs, stats = c("elpd", "auc"))$selection
  df["rho"] <- rho
  df["n"] <- N
  df["loo_wts"] <- exp(df["elpd"]) / sum(exp(df["elpd"]))

  # write the table
  setwd("/scratch/work/mclatcy1/projpred-workflow")
  csv_name <- paste0("./data/binomial/",N,"n,",rho,"rho.csv")
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
  N = options[,2], 
  mc.cores = Sys.getenv('SLURM_CPUS_PER_TASK')
)

# experiment_old <- function (N, rho) {
#   # simulate data
#   data <- .generate(N, rho)
#   df <- as.data.frame(data)
  
#   # fit binomial family GLM
#   fit <- brms::brm(
#     formula = "y ~ .", 
#     family = bernoulli(link = "logit"),
#     prior = brms::set_prior("horseshoe(3)"), 
#     data = df
#   )

#   # perform projection predictive inference
#   cv <- cv_varsel(
#     fit, 
#     method = "forward", 
#     validate_search = TRUE, 
#     nterms_max = 25
#   )

#   # compute reference model criteria on independent test data
#   newdata <- .generate(N, rho)
#   new_llk <- log_lik(fit, newdata = newdata)
#   ref.preds <- predict(fit, newdata = newdata)[, "Estimate"]
#   ref.elpd <- loo::elpd(new_llk)$estimates["elpd", "Estimate"]
#   ref.acc <- .acc(y_true = indep_data$y, y_pred = ref.preds)
#   ref.auc <- .auc(y_true = indep_data$y, y_pred = ref.preds)
  
#   # initialise summary dataframe
#   summary <- data.frame(matrix(nrow = 0, ncol = 11))

#   # iterate over solution path
#   for (size in 0:length(cv$solution_terms)) {

#     # generate new data
#     newdata <- .generate(N, rho)

#     # project onto subset
#     sel_terms <- head(cv$search_path$solution_terms, size)
#     prj <- project(
#       fit,
#       solution_terms = sel_terms
#     )
    
#     # compute predictions on independent test data
#     prj_preds <- proj_linpred(prj, newdata = newdata, integrated = TRUE)
    
#     # compute submodel criteria
#     acc <- .acc(y_true = newdata$y, y_pred = prj_preds$pred)
#     acc.diff <- ref.acc - acc
#     auc <- .auc(y_true = newdata$y, y_pred = prj_preds$pred)
#     auc.diff <- ref.auc - auc
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
#       acc = rmse,
#       acc.diff = acc.diff,
#       auc = rmse,
#       auc.diff = auc.diff
#     )
#     summary <- rbind(summary, submdl_summary)
#   }

#   # compute LOO weights across submodels
#   summary["loo_wts"] <- exp(summary["elpd"]) / sum(exp(summary["elpd"]))

#   # write the table
#   csv_name <- paste0("./data/binomial/",N,"n,",rho,"rho.csv")
#   write.table(summary, file = csv_name, sep = ",", row.names = FALSE)
# }
