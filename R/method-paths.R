# load necessary libraries
library(tidyverse)
library(brms)
library(projpred)
library(parallel)

# load data generating scripts
source("./R/generate.R")

# define experiment
experiment <- function(rho, method){
  print(rho)
  print(method)
  # define size of data to use in experiments
  n <- 200
  
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
  vs <- varsel(
    fit, 
    method = method, 
    nterms_max = 100
  )
  
  paths <- vs$solution_terms
  
  for (i in 1:10) {
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
    vs <- varsel(
      fit, 
      method = method, 
      nterms_max = 100
    )
    
    # append selection paths
    paths <- cbind(paths, vs$solution_terms)
  }
  
  # helper function to compute occurrences of covariate rankings
  .count_occs <- function(x) {
    occs <- which(paths == x, arr.ind = TRUE)
    while (nrow(occs) < ncol(paths)) {
      pad <- data.frame(row = 0, col = 0)
      occs <- rbind(occs, pad)
    }
    occs[occs == 0] <- NA
    res <- as.data.frame(table(occs[,"row"]))
    print(res)
    res$Freq <- res$Freq / sum(res$Freq)
    return (res)
  }
  
  els <- unique(as.vector(as.matrix(paths)))
  occs <- purrr::map(els, .count_occs)
  names(occs) <- els
  res <- do.call("rbind", occs)
  
  res_names <- sub("(.*?..*?).(.*?)", "\\1", rownames(res))
  res_names <- sub("\\..*", "", res_names)
  res["name"] <- res_names
  res["rho"] <- rho
  res["method"] <- method
  res <- rename(res, rank=Var1, freq=Freq)
  return(res)
}

# vary both the number of data observations and the correlation between parameters
rhos <- c(0, 0.5, 0.9)
methods <- c("forward", "L1")
options <- expand.grid(rho=rhos, method=methods)

# perform experiment
res <- parallel::mcMap(
  f = experiment, 
  rho = options[,1], 
  method = options[,2], 
  mc.cores = Sys.getenv('SLURM_CPUS_PER_TASK')
)

# concatenate experiment results
df <- do.call("rbind", res)

# write the table
csv_name <- paste0("./data/paths/method/method_paths.csv")
ff <- file(csv_name, open="w")
write.table(df, file = ff, sep = ",", row.names = FALSE)
close(ff)
