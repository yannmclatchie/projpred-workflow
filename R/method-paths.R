# load necessary libraries
library(brms)
library(projpred)
library(parallel)
library(ggplot2)
library(geomtextpath)
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
ref_loo <- loo(fit)
refm_elpd <- ref_loo$estimates["elpd_loo", "Estimate"]

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

# plot the solution paths of both methods
vs_forward_no_cv <- cv_varsel(
  fit,
  method = "forward",
  validate_search = FALSE,
  nclusters_pred = 10,
  seed = SEED
)
vs_l1_no_cv <- cv_varsel(
  fit,
  method = "l1",
  validate_search = FALSE,
  nclusters_pred = 10,
  seed = SEED
)
forward_df <- vs_forward_no_cv$summary %>%
  select(c(size, elpd.loo, se)) %>%
  mutate(method = "Forward")
l1_df <- vs_l1_no_cv$summary %>%
  select(c(size, elpd.loo, se)) %>%
  mutate(method = "L1")

( gg_path <- forward_df %>%
  rbind(l1_df) %>%
  #filter(size > 0) %>%
  ggplot(
    aes(
      x = size, 
      y = elpd.loo,
      ymin = elpd.loo - se, 
      ymax = elpd.loo + se,
      colour = method,
      label = method
    ) 
  ) + 
    geom_hline(yintercept = refm_elpd, colour = "red", linetype = "longdash") + 
    geom_pointrange(position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0)) +
    geom_labelpath(
      aes(vjust = method, hjust = method),
      position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0),
      text_smoothing = 75
    ) +
    annotate("text", x = 10, y = refm_elpd - 5, colour = "red", label = "Reference model elpd") +
    ylab("$elpd$") +
    xlab("Model size") +
    scale_color_manual(values = c("black", "blue")) +
    scale_vjust_manual(values = c(-3, 3)) + 
    scale_hjust_manual(values = c(0.3, 0.3)) +
    scale_y_continuous(limits = c(-760, -700)) +
    theme_bw() + 
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position="none") )

source("./R/aux/aux_plotting.R")
save_tikz_plot(plot = gg_path,
               filename = "./tex/forward_l1_paths.tex",
               width = 3.5)

# plot the over-fitting search path
forward_df_cv <- vs_forward$summary %>%
  select(c(size, elpd.loo, se)) %>%
  mutate(cv_search = "cross-validated search")
( p_elpd <- forward_df %>%
  mutate(cv_search = "full-data search") %>%
  select(-method) %>%
  rbind(forward_df_cv) %>%
  ggplot(
    aes(
      x = size, 
      y = elpd.loo,
      ymin = elpd.loo - se, 
      ymax = elpd.loo + se,
      colour = cv_search,
      label = cv_search
    )
  ) +
    geom_hline(yintercept = refm_elpd, colour = "red", linetype = "longdash") + 
    geom_pointrange(position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0)) +
    geom_labelpath(
      aes(vjust = cv_search, hjust = cv_search),
      position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0),
      text_smoothing = 75
    ) +
    annotate("text", x = 0.5, y = refm_elpd + 5, colour = "red", label = "Reference model elpd") +
    ylab("$elpd$") +
    xlab("Model size") +
    scale_color_manual(values = c("black", "blue")) +
    scale_vjust_manual(values = c(5, -5)) + 
    scale_hjust_manual(values = c(0.3, 0.5)) +
    scale_y_continuous(limits = c(-760, -700)) +
    theme_bw() + 
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position="none") )
save_tikz_plot(plot = p_elpd,
               filename = "./tex/cv_forward_paths.tex",
               width = 3.5)

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

save_tikz_plot(plot = gg_fw_l1_fixedY,
               filename = "./tex/pct_solution_terms_cv_forward_l1_fixedY.tex",
               width = 7)
