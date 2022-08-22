library("bayesplot")
library("ggplot2")
library("brms")

# load data generating scripts
source("./R/generate.R")

data <- .generate(0.5, 100)
fit <- brms::brm(
  formula = "y ~ .", 
  family = "gaussian", 
  prior = brms::set_prior("horseshoe(11)"), 
  data = data
)
prior <- get_prior(formula = "y ~ .", 
                   family = "gaussian", 
                   prior = brms::set_prior("horseshoe(11)"), 
                   data = data)
unique(prior[, "class"])
prior[prior[, "class"] == "sigma", ]

bayesplot::mcmc_pairs(fit)
rhats <- bayesplot::rhat(fit)
bayesplot::color_scheme_set("brightblue")
p_rhat <- bayesplot::mcmc_rhat(rhats)
p_rhat <- p_rhat + xaxis_title(size = 30) + xaxis_text(size = 30) + legend_text(size = 30)


bayesplot::color_scheme_set("brightblue")
yrep <- brms::posterior_predict(fit, draws = 500)
p_ppc_dens <- bayesplot::ppc_dens_overlay(data$y, yrep[1:50, ])
p_ppc_dens <- p_ppc_dens + xaxis_title(size = 30) + xaxis_text(size = 30) + legend_text(size = 30)
p_ppc_hist <- bayesplot::ppc_hist(data$y, yrep[1:5, ])
p_ppc_hist <- p_ppc_hist + xaxis_title(size = 30) + xaxis_text(size = 30) + legend_text(size = 30)
p_ppc_mean <- bayesplot::ppc_stat(data$y, yrep, stat = "mean")
p_ppc_mean <- p_ppc_mean + xaxis_title(size = 30) + xaxis_text(size = 30) + legend_text(size = 30)
p_ppc_mean

# save plot
ggsave(plot = p_rhat, filename = "./img/ppc_rhat.pdf")
ggsave(plot = p_ppc_dens, filename = "./img/ppc_dens.pdf")
ggsave(plot = p_ppc_hist, filename = "./img/ppc_hist.pdf")
ggsave(plot = p_ppc_mean, filename = "./img/ppc_mean.pdf")
