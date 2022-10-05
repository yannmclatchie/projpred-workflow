# load necessary libraries
library(tidyverse)
library(brms)
library(projpred)
library(parallel)
library(scam)

# load data generating scripts
source("./R/generate.R")

num_runs <- 2
data <- experiment(0.5, 200)
data

means = data %>% 
  group_by(method, size, n, rho) %>% 
  summarize(
    mean_diff = mean(elpd.diff)
  )

data %>%
  ggplot() +
  geom_line(
    data = filter(data, method == "raw"), 
    aes(x = size, y = elpd.diff, group = run, colour = "Raw"), 
    alpha = 0.2
  ) +
  geom_line(
    data = filter(data, method == "smoothed"), 
    aes(x = size, y = elpd.diff, group = run, colour = "Smoothed"), 
    alpha = 0.2
  ) +
  geom_line(
    data = filter(means, method == "raw"),
    aes(x = size, y = mean_diff, color = 'Raw'),
    size = 0.9
  ) +
  geom_line(
    data = filter(means, method == "smoothed"),
    aes(x = size, y = mean_diff, color = 'Smoothed'),
    size = 0.9
  ) +
  scale_colour_manual(
    "", 
    breaks = c("Raw", "Smoothed"),
    values = c("#EE6677", "#4477AA")) +
  facet_grid(rows = vars(n), cols = vars(rho)) + 
  ylab('ELPD difference')+
  xlab('Model size')+
  theme_classic()


geom_line(data = filter(means, name == strategy),
          aes(x = size, y = mean_mlpd_loo, color = 'Train (MLPD LOO)'),
          size=0.9) +
  geom_line(data = filter(means, name == strategy),
            aes(x = size, y = mean_mlpd_test, color = 'Test (MLPD)'),
            size=0.9) +
  scale_colour_manual("", 
                      breaks = c("Train (MLPD LOO)", "Test (MLPD)"),
                      values = c("red", "black")) +
  facet_grid(rows = vars(n), cols = vars(rho))

# define experiment
experiment <- function(rho, n){
  
  df <- data.frame(matrix(ncol = 13, nrow = 0))
  
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
      nterms_max = 25
    )
    
    # compute submodel summaries
    run_df <- summary(vs, stats = c("mlpd", "elpd"))$selection
    run_df["rho"] <- rho
    run_df["n"] <- n
    run_df["run"] <- run
    run_df["method"] <- "raw"
    run_df
    
    # build spline dataframe
    spline_df <- run_df
    spline_df["method"] <- "smoothed"
    
    # fit different splines to the projpred results
    spline <- scam(elpd.diff/elpd.diff.se ~ s(size, bs = "mpi"), data = spline_df)
    
    # update dataframe
    spline_df <- spline_df %>%
      mutate(elpd.diff = spline$fit*spline_df[,'elpd.diff.se'],
             elpd.diff.se = sqrt(spline$sig2)*spline_df[,'elpd.diff.se'])
    run_df <- rbind(run_df, spline_df)
    df <- rbind(df, run_df)
  }
  
  return (df)
}

summary <- experiment(rho = 0.5, n = 200)
summary

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