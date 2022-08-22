library(projpred)
library(tidyverse)
library(brms)
library(parallel)

# define experiment
experiment <- function(rho, n){
  
  # read data
  csv_name <- paste0("./data/generated/",n,"n,",rho,"rho.csv")
  data <- vroom::vroom(csv_name)
  
  # define number of parameters in reference model
  p <- 100

  # fit the BRMS model
  fit <- brm(
    formula = "y ~ .", 
    family = "gaussian", 
    prior = set_prior("horseshoe(11)"), 
    data = data
  )
  
  # using 50 clusters
  ref_model_50_clusters <- cv_varsel(
    fit, 
    nterms_max = p,
    ndraws = NULL,
    nclusters = 50,
    ndraws_pred = 400,
    nclusters_pred = NULL,
    validate_search = FALSE
  )
  summary <- ref_model_50_clusters$summary
  summary["nclusters"] = 50
  summary["rho"] = rho
  summary["n"] = n
  
  # using recommended parameters
  ref_model_20_clusters <- cv_varsel(
    fit, 
    nterms_max = p,
    ndraws = NULL,
    nclusters = 20,
    ndraws_pred = 400,
    nclusters_pred = NULL,
    validate_search = FALSE
  )
  summary20 <- ref_model_20_clusters$summary
  summary20["nclusters"] = 20
  summary20["rho"] = rho
  summary20["n"] = n
  summary <- rbind(summary, summary20)
  
  # using 10 clusters
  ref_model_10_clusters <- cv_varsel(
    fit, 
    nterms_max = p,
    ndraws = NULL,
    nclusters = 10,
    ndraws_pred = 400,
    nclusters_pred = NULL,
    validate_search = FALSE
  )
  summary10 <- ref_model_10_clusters$summary
  summary10["nclusters"] = 10
  summary10["rho"] = rho
  summary10["n"] = n
  summary <- rbind(summary, summary10)
  
  # using 5 clusters
  ref_model_5_clusters <- cv_varsel(
    fit, 
    nterms_max = p,
    ndraws = NULL,
    nclusters = 5,
    ndraws_pred = 400,
    nclusters_pred = NULL,
    validate_search = FALSE
  )
  summary5 <- ref_model_5_clusters$summary
  summary5["nclusters"] = 5
  summary5["rho"] = rho
  summary5["n"] = n
  summary <- rbind(summary, summary5)

  # append the csv
  write.table(
    summary, 
    file = "./data/clusters.csv", 
    append = TRUE, 
    sep = ",",
    row.names = FALSE,
    col.names = FALSE
  )
}

# initialise the data
summary <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(summary) <- c(
    "size", 
    "solution_terms", 
    "elpd.loo", 
    "se", 
    "diff", 
    "diff.se", 
    "nclusters", 
    "rho", 
    "n"
)
write.table(summary, file = "./data/clusters.csv", sep = ",")

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
