# load necessary libraries
library(dplyr)
library(brms)
library(projpred)
library(ggplot2)

# get data
data("df_gaussian", package = "projpred")
dat_gauss <- data.frame(y = df_gaussian$y, df_gaussian$x)

# define reference model
( D <- sum(grepl("^X", names(dat_gauss))) )
# guess of the number of relevant
p0 <- 5
# number of observations:
N <- nrow(dat_gauss)
# hyperprior scale for tau
tau0 <- p0 / (D - p0) * 1 / sqrt(N)

# fit the reference model
refm_fit <- brm(
  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 +
    X15 + X16 + X17 + X18 + X19 + X20,
  family = gaussian(),
  data = dat_gauss,
  prior = set_prior(horseshoe(scale_global = tau0)),
  seed = 2052109, refresh = 0
)

# run projpred
cvvs <- cv_varsel(
  refm_fit,
  validate_search = FALSE,
  nclusters_pred = 20,
  nterms_max = D,
  seed = 411183
)

# repeat the search path up to the point of most bulge
vs <- cv_varsel(
  refm_fit,
  validate_search = TRUE,
  nterms_max = 7,
  seed = 411183
)

# compute submodel summaries
sel_df <- summary(cvvs, stats = c("elpd"))$selection
sel_df[, "cv_search"] = FALSE

# append validated search results to the data
cv_sel_df <- summary(vs, stats = c("elpd"))$selection
cv_sel_df[, "cv_search"] = TRUE
sel_df <- rbind(sel_df, cv_sel_df)

# compute the reference model's LOO-CV predictive performance
refm_loo <- loo(refm_fit)
refm_elpd <- refm_loo$estimates["elpd_loo", "Estimate"]

# plot the over-fitting search path
p_elpd <- sel_df %>%
  ggplot(
    aes(
      x = size, 
      y = elpd.loo,
      ymin = elpd.loo - se, 
      ymax = elpd.loo + se,
      colour = cv_search
    )
  ) +
  geom_hline(yintercept = refm_elpd, colour = "red", linetype = "longdash") + 
  geom_vline(xintercept = 7, colour = "grey", linetype = "longdash") + 
  geom_pointrange(position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0)) +
  geom_line(position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0)) +
  annotate("text", x = 2, y = refm_elpd + 5, colour = "red", label = "Reference model elpd") +
  ylab("$elpd$") +
  xlab("Model size") +
  scale_color_manual(values = c("black", "blue")) +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position="none")
p_elpd

# save plot to tikz
source("./R/aux/aux_plotting.R")
save_tikz_plot(p_elpd, width = 6, filename = "./tex/elpd-bulge.tex")

