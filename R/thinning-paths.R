# load necessary libraries
library(tidyverse)
library(brms)
library(projpred)
library(parallel)

# load data generating scripts
source("./R/generate.R")

# define experiment
experiment <- function(rho){
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
  vs_draws <- varsel(
    fit, 
    nclusters = NULL,
    ndraws = 50,
    method = "forward", 
    nterms_max = 100
  )
  vs_clusters <- varsel(
    fit, 
    nclusters = 20,
    ndraws = NULL,
    method = "forward", 
    nterms_max = 100
  )
  
  draws_paths <- vs_draws$solution_terms
  clusters_paths <- vs_clusters$solution_terms
  
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
    vs_draws <- varsel(
      fit, 
      nclusters = NULL,
      ndraws = 50,
      method = "forward", 
      nterms_max = 100
    )
    vs_clusters <- varsel(
      fit, 
      nclusters = 20,
      ndraws = NULL,
      method = "forward", 
      nterms_max = 100
    )
    
    draws_paths <- cbind(draws_paths, vs_draws$solution_terms)
    clusters_paths <- cbind(clusters_paths, vs_clusters$solution_terms)
  }
  
  # helper function to compute occurrences of covariate rankings
  .count_draw_occs <- function(x) {
    occs <- which(draws_paths == x, arr.ind = TRUE)
    while (nrow(occs) < ncol(draws_paths)) {
      pad <- data.frame(row = 0, col = 0)
      occs <- rbind(occs, pad)
    }
    occs[occs == 0] <- NA
    res <- as.data.frame(table(occs[,"row"]))
    print(res)
    res$Freq <- res$Freq / sum(res$Freq)
    return (res)
  }
  .count_cluster_occs <- function(x) {
    occs <- which(clusters_paths == x, arr.ind = TRUE)
    while (nrow(occs) < ncol(clusters_paths)) {
      pad <- data.frame(row = 0, col = 0)
      occs <- rbind(occs, pad)
    }
    occs[occs == 0] <- NA
    res <- as.data.frame(table(occs[,"row"]))
    print(res)
    res$Freq <- res$Freq / sum(res$Freq)
    return (res)
  }
  
  draw_els <- unique(as.vector(as.matrix(draws_paths)))
  cluster_els <- unique(as.vector(as.matrix(clusters_paths)))

  draw_occs <- purrr::map(draw_els, .count_draw_occs)
  cluster_occs <- purrr::map(cluster_els, .count_cluster_occs)

  names(draw_occs) <- draw_els
  names(cluster_occs) <- cluster_els

  draw_res <- do.call("rbind", draw_occs)
  cluster_res <- do.call("rbind", cluster_occs)
  
  draw_res_names <- sub("(.*?..*?).(.*?)", "\\1", rownames(draw_res))
  draw_res_names <- sub("\\..*", "", draw_res_names)
  draw_res["name"] <- draw_res_names
  draw_res["rho"] <- rho
  draw_res["thinning"] <- "draws"
  draw_res <- rename(draw_res, rank=Var1, freq=Freq)

  cluster_res_names <- sub("(.*?..*?).(.*?)", "\\1", rownames(cluster_res))
  cluster_res_names <- sub("\\..*", "", cluster_res_names)
  cluster_res["name"] <- cluster_res_names
  cluster_res["rho"] <- rho
  cluster_res["thinning"] <- "clusters"
  cluster_res <- rename(cluster_res, rank=Var1, freq=Freq)

  res <- rbind(draw_res, cluster_res)
  return(res)
}

# vary both the number of data observations and the correlation between parameters
rhos <- c(0, 0.5, 0.9)

# perform experiment
res <- parallel::mcMap(
  f = experiment, 
  rho = rhos,
  mc.cores = Sys.getenv('SLURM_CPUS_PER_TASK')
)

# concatenate experiment results
df <- do.call("rbind", res)

# write the table
setwd("/scratch/work/mclatcy1/projpred-workflow")
csv_name <- paste0("./data/paths/thinning/thinning_paths.csv")
ff <- file(csv_name, open="w")
write.table(df, file = ff, sep = ",", row.names = FALSE)
close(ff)