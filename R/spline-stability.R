# check if scam if installed
if (!require("scam")) install.packages("scam", repos = "http://cran.us.r-project.org"); library(scam)

# load necessary libraries
library(dplyr)
library(brms)
library(projpred)
library(parallel)
library(gamm4)

# load data generating scripts
source("./R/generate.R")

num_runs <- 10

# define experiment
experiment <- function(rho, n){
  
  df <- data.frame(matrix(ncol = 15, nrow = 0))
  
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
    
    # # select model size for each method
    # loo.elpd.delta.size <- run_df %>%
    #   filter( (diff >= -4) ) %>%
    #   arrange(desc(as.numeric(size))) %>%
    #   slice(n = n()) %>%
    #   dplyr::select(size)
    # spline.elpd.delta.size <- run_df %>%
    #   filter( (spline.elpd.diff >= -4) ) %>%
    #   arrange(desc(as.numeric(size))) %>%
    #   slice(n = n()) %>%
    #   dplyr::select(size)
    # mono.spline.elpd.delta.size <- run_df %>%
    #   filter( (mono.spline.elpd.diff >= -4) ) %>%
    #   arrange(desc(as.numeric(size))) %>%
    #   slice(n = n()) %>%
    #   dplyr::select(size)
    # loo.elpd.ci.size <- run_df %>%
    #   filter( (elpd + se >= ref.loo.elpd) ) %>%
    #   arrange(desc(as.numeric(size))) %>%
    #   slice(n = n()) %>%
    #   dplyr::select(size)
    # spline.elpd.ci.size <- run_df %>%
    #   filter( (spline.elpd + spline.elpd.se >= ref.loo.elpd) ) %>%
    #   arrange(desc(as.numeric(size))) %>%
    #   slice(n = n()) %>%
    #   dplyr::select(size)
    # mono.spline.elpd.ci.size <- run_df %>%
    #   filter( (mono.spline.elpd + mono.spline.elpd.se >= ref.loo.elpd) ) %>%
    #   arrange(desc(as.numeric(size))) %>%
    #   slice(n = n()) %>%
    #   dplyr::select(size)
    # loo.elpd.switch.size <- min(loo.elpd.ci.size, loo.elpd.delta.size)
    # spline.elpd.switch.size <- min(spline.elpd.ci.size, spline.elpd.delta.size)
    # mono.spline.elpd.switch.size <- min(
    #   mono.spline.elpd.ci.size, mono.spline.elpd.delta.size
    # )
    # sel_df <- data.frame(
    #   loo.elpd.delta = loo.elpd.delta.size,
    #   spline.elpd.delta = spline.elpd.delta.size,
    #   mono.spline.elpd.delta = mono.spline.elpd.delta.size,
    #   loo.elpd.ci = loo.elpd.ci.size,
    #   spline.elpd.ci = spline.elpd.ci.size,
    #   mono.spline.elpd.ci = mono.spline.elpd.ci.size,
    #   loo.elpd.switch = loo.elpd.switch.size,
    #   spline.elpd.switch = spline.elpd.switch.size,
    #   mono.spline.elpd.switch = mono.spline.elpd.switch.size
    # )
    
    # # melt data and add run details
    # res <- reshape2::melt(sel_df, variable.name = "procedure", value.name = "size")
    # res["n"] <- n
    # res["rho"] <- rho
    
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
