library(tidyverse)
library(brms)
library(projpred)
library(tictoc)

# load data generating scripts
source("./R/generate.R")

time_projpred <- function (method, validate_search, nterms_max, nclusters) {
  message(method, validate_search, nterms_max, nclusters)

  # initialise dataframe
  df <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(df) <- c("method", "validate_search", "nterms_max", "nclusters", "time")
  
  # run experiment 5 times
  for (i in 1:5) {
    # simulate data
    data <- .generate(0.5, 200)

    # fit the BRMS model
    fit <- brm(
      formula = "y ~ .", 
      family = "gaussian", 
      prior = set_prior("horseshoe(11)"), 
      data = data
    )

    # start timer
    tic()
    
    # run projpred
    cv_varsel(
      fit, 
      method = method, 
      validate_search = validate_search, 
      nterms_max = nterms_max,
      nclusters = nclusters
    )
    
    # log time taken
    exectime <- toc()
    exectime <- exectime$toc - exectime$tic
    
    # append to dataframe
    row <- data.frame(
      method = method, 
      validate_search = validate_search, 
      nterms_max = nterms_max,
      nclusters = nclusters,
      time = exectime
      )
    message(row)
    df <- rbind(df, row)
  }
  return (df)
}

# vary the timing options
methods <- c("forward", "L1")
validations<- c(TRUE, FALSE)
nterms <- seq(10, 100, length.out = 3)
nclusters <- seq(20, 2000, length.out = 3)
options <- rbind(
  reshape::expand.grid.df(
    data.frame(
      validations=TRUE, 
      nterms=100, 
      nclusters=20
    ),
    data.frame(methods=methods)
  ),
  reshape::expand.grid.df(
    data.frame(
      validations=TRUE,
      methods="forward",
      nclusters=20
    ),
    data.frame(nterms=nterms)
  ),
  reshape::expand.grid.df(
    data.frame(
      validations=TRUE, 
      nterms=100, 
      methods="forward"
    ),
    data.frame(nclusters=nclusters)
  ),
  reshape::expand.grid.df(
    data.frame(
      methods="forward",
      nterms=100, 
      nclusters=20
    ),
    data.frame(validations=validations)
  )
) %>%
  distinct(.keep_all = TRUE) # keep only distinct combinations
message(nrow(options))

# perform experiment
res <- parallel::mcMap(
  f = time_projpred, 
  method = options[,"methods"],
  validate_search = options[,"validations"], 
  nterms_max = options[,"nterms"], 
  nclusters = options[,"nclusters"],
  mc.cores = Sys.getenv('SLURM_CPUS_PER_TASK')
)

# concatenate experiment results
df <- do.call("rbind", res)

# save the concatenated results
setwd("/scratch/work/mclatcy1/projpred-workflow")
csv_name <- paste0("./data/timing.csv")
ff <- file(csv_name, open="w")
write.table(df, file = ff, sep = ",", row.names = FALSE)
close(ff)
