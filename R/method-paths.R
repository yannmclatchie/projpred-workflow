# load necessary libraries
library(brms)
library(projpred)
library(parallel)
library(ggplot2)
library(patchwork)

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
#vs_forward <- cv_varsel(
#  fit,
#  method = "forward",
#  nclusters_pred = 10,
#  seed = SEED
#)
#saveRDS(vs_forward, "data/vs_forward.rds")
vs_forward <- readRDS("data/vs_forward.rds")
#vs_l1 <- cv_varsel(
#  fit,
#  method = "l1",
#  nclusters_pred = 10,
#  seed = SEED
#)
#saveRDS(vs_l1, "data/vs_l1.rds")
vs_l1 <- readRDS("data/vs_l1.rds")

# plot the stability of selection process
source("./R/aux/projpredpct.R")
source("./R/aux/gg_pct_solution_terms_cv.R")
( gg_fw <- gg_pct_solution_terms_cv(vs_forward) 
  & labs(subtitle = "Forward search"))
( gg_l1 <- gg_pct_solution_terms_cv(vs_l1) 
  & theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) 
  & labs(subtitle = "$L_1$ search"))

( gg_fw_l1 <- gg_fw | gg_l1 )
y_fw_order <- ggplot_build(gg_fw)$layout$panel_scales_y[[1]]$range$range
gg_fw_l1_fixedY <- (gg_fw + gg_l1) & scale_y_discrete(limits = y_fw_order)
gg_fw_l1_fixedY

source("./R/aux/aux_plotting.R")
save_tikz_plot(plot = gg_fw_l1_fixedY,
               filename = "./tex/pct_solution_terms_cv_forward_l1_fixedY.tex",
               width = 7)
