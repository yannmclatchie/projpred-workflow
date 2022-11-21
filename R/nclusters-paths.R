# load necessary libraries
library(brms)
library(projpred)
library(parallel)

# set seed
SEED <- 1513306866

# load data
df <- read.table("./data/bodyfat/bodyfat.txt", header = T, sep = ";")
df[,4:19] <- scale(df[,4:19])
# no-one can have 0% body fat
df <- df[df$siri>0,]
df <- as.data.frame(df)
( n <- nrow(df) )

# define the covariates and variate
pred <- c("age", "weight", "height", "neck", "chest", "abdomen", "hip",
          "thigh", "knee", "ankle", "biceps", "forearm", "wrist")
target <- "siri"
formula <- formula(paste("siri ~", paste(pred, collapse = " + ")))

# fit the reference model
options(mc.cores = detectCores(logical = FALSE))
r2d2_prior <- set_prior(R2D2(mean_R2 = 0.3, prec_R2 = 3))
# indep_gauss_prior <- set_prior("normal(0, 2.5)", class = "b")
fit <- brm(
  formula,
  data = df,
  prior = r2d2_prior, # prior = indep_gauss_prior,
  seed = SEED,
  refresh = 0,
  file = "bodyfat_fit",
  file_refit = "on_change"
)

# perform projpred with forward and L1 search
vs_10_clust <- cv_varsel(
  fit,
  method = "forward",
  nclusters = 10,
  nclusters_pred = 10,
  seed = SEED
)
saveRDS(vs_10_clust, "vs_10_clust.rds")
# vs_10_clust <- readRDS("vs_10_clust.rds")
vs_50_clust <- cv_varsel(
  fit,
  method = "forward",
  nclusters = 50,
  nclusters_pred = 10,
  seed = SEED
)
saveRDS(vs_50_clust, "vs_50_clust.rds")
# vs_50_clust <- readRDS("vs_50_clust.rds")
vs_200_clust <- cv_varsel(
  fit,
  method = "forward",
  nclusters = 200,
  nclusters_pred = 10,
  seed = SEED
)
saveRDS(vs_200_clust, "vs_200_clust.rds")
# vs_200_clust <- readRDS("vs_200_clust.rds")

# plot the stability of selection process
source("./R/aux/projpredpct.R")
source("./R/aux/gg_pct_solution_terms_cv.R")
( gg_10 <- gg_pct_solution_terms_cv(vs_10_clust) )
( gg_50 <- gg_pct_solution_terms_cv(vs_50_clust) )
( gg_200 <- gg_pct_solution_terms_cv(vs_200_clust) )
library(patchwork)
( gg_all <- gg_10 / gg_50 / gg_200 )
y_10_order <- ggplot_build(gg_10)$layout$panel_scales_y[[1]]$range$range
( gg_all_fixedY <- (gg_10 / gg_50 / gg_200) & scale_y_discrete(limits = y_10_order) )
source("./R/aux/aux_plotting.R")
save_tikz_plot(plot = gg_10,
               filename = "./tex/pct_solution_terms_cv_10.tex",
               width = 6)
save_tikz_plot(plot = gg_50,
               filename = "./tex/pct_solution_terms_cv_50.tex",
               width = 6)
save_tikz_plot(plot = gg_200,
               filename = "./tex/pct_solution_terms_cv_200.tex",
               width = 6)
save_tikz_plot(plot = gg_all,
               filename = "./tex/pct_solution_terms_cv_10_50_200.tex",
               width = 6)
save_tikz_plot(plot = gg_all_fixedY,
               filename = "./tex/pct_solution_terms_cv_all_fixedY.tex",
               width = 6)
