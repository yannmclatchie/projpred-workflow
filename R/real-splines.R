# check if scam if installed
if (!require("scam")) install.packages("scam"); library("scam")

# load necessary libraries
library(tidyverse)
library(brms)
library(projpred)
library(parallel)
library(gamm4)
library(R.matlab)

experiment <- function(data) {
  if (data == "GLIOMA") {
    root <- "./data/benchmarks/GLIOMA.mat"
    family <- "categorical"
  } else if (data == "leukemia") {
    root <- "./data/benchmarks/leukemia.mat"
    family <- "bernoulli"
  } else if (data == "colon") {
    root <- "./data/benchmarks/colon.mat"
    family <- "bernoulli"
  } else if (data == "prostate") {
    root <- "./data/benchmarks/Prostate_GE.mat"
    family <- "bernoulli"
  } else {
    warning("Unavailable data")
  }
  
  # read matrix
  mat <- readMat(root)
  
  # build dataframe
  data <- data.frame(y = mat$Y)
  num_features <- length(mat$X[1,])
  for (i in 1:num_features) {
    Newcolname <- paste0("X.",i) 
    data[Newcolname] <- mat$X[,i]
  }
  
  # fit the reference model
  fit <- brm(
    formula = "y ~ .", 
    family = family, 
    prior = set_prior(R2D2(mean_R2 = 0.5, prec_R2 = 50, cons_D2 = 1, autoscale = TRUE)), 
    data = data
  )
  
  # perform projections
  vs <- varsel(
    fit, 
    method = "forward", 
    nterms_max = 25
  )
  
  # compute submodel summaries
  sel_df <- summary(vs, stats = c("elpd"))$selection
  
  # extract reference model ELPD
  ref.loo <- loo(fit)
  sel_df["ref.loo.elpd"] <- ref.loo$estimates["elpd_loo", "Estimate"]
  
  # fit different splines to the projpred results
  spline <- gam(elpd/se ~ s(size), data = sel_df)
  mono.spline <- scam(elpd/se ~ s(size, bs = "mpi"), data = sel_df)
  
  # update dataframe
  sel_df <- sel_df %>%
    mutate(
      spline.elpd = spline$fit*run_df[,'se'],
      spline.elpd.se = sqrt(spline$sig2)*run_df[,'se'],
      mono.spline.elpd = mono.spline$fit*run_df[,'se'],
      mono.spline.elpd.se = sqrt(mono.spline$sig2)*run_df[,'se']
    )
  sel_df["data"] <- data
  return(sel_df)
}
