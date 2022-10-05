# install.packages("remotes")
# remotes::install_github("n-kall/priorsense")

library("bayesplot")
library("ggplot2")
library("brms")
library("priorsense")

# load data generating scripts
source("./R/generate.R")

data <- .generate(0.5, 100)
names(data) <- gsub("\\.", "", names(data))
fit <- brms::brm(
  formula = "y ~ .", 
  family = "gaussian", 
  prior = brms::set_prior("horseshoe(11)"), 
  data = data
)
priorsense::powerscale_sensitivity(fit)

pss <- powerscale_sequence(fit)
p_pss <- powerscale_plot_ecdf(pss, variables = c("b_Intercept", "b_X1", "b_X2", "b_X3", "b_X4", "b_X5")) + theme_classic() + theme(title = element_blank())
p_pss

prior <- get_prior(formula = "y ~ .", 
                   family = "gaussian", 
                   prior = brms::set_prior("horseshoe(11)"), 
                   data = data)
unique(prior[, "class"])
prior[prior[, "class"] == "sigma", ]


loo <- loo(fit, save_psis = TRUE, cores = 2)
plot(loo, label_points = TRUE)

ratios <- neff_ratio(fit)
p_neff <- mcmc_neff(ratios, size = 2) + xaxis_title(size = 30) + xaxis_text(size = 30) + legend_text(size = 30)
p_neff

bayesplot::mcmc_pairs(fit)
rhats <- bayesplot::rhat(fit)
bayesplot::color_scheme_set("brightblue")
p_rhat <- bayesplot::mcmc_rhat(rhats)
p_rhat <- p_rhat + xaxis_title(size = 30) + xaxis_text(size = 30) + legend_text(size = 30)


yrep <- brms::posterior_predict(fit, draws = 500)
p_ppc_dens <- bayesplot::ppc_dens_overlay(data$y, yrep[1:50, ])
p_ppc_dens <- p_ppc_dens + xaxis_title(size = 30) + xaxis_text(size = 30) + legend_text(size = 30)
p_ppc_hist <- bayesplot::ppc_hist(data$y, yrep[1:5, ])
p_ppc_hist <- p_ppc_hist + xaxis_title(size = 30) + xaxis_text(size = 30) + legend_text(size = 30)
p_ppc_mean <- bayesplot::ppc_stat(data$y, yrep, stat = "mean")
p_ppc_mean <- p_ppc_mean + xaxis_title(size = 30) + xaxis_text(size = 30) + legend_text(size = 30)
p_ppc_mean

ppc_dens(data$y, yrep[100:102, ])
p_ecdf <- ppc_ecdf_overlay(data$y, yrep[1:50, ]) + xaxis_text(size = 30) + yaxis_text(size = 30) + legend_text(size = 30)
p_ecdf

# save plot
ggsave(plot = p_pss, filename = "./img/pss.pdf", width = 16, height = 9, units = "cm")
ggsave(plot = p_neff, filename = "./img/ppc_neff.pdf")
ggsave(plot = p_rhat, filename = "./img/ppc_rhat.pdf")
ggsave(plot = p_ppc_dens, filename = "./img/ppc_dens.pdf")
ggsave(plot = p_ecdf, filename = "./img/ppc_ecdf.pdf")
ggsave(plot = p_ppc_hist, filename = "./img/ppc_hist.pdf")
ggsave(plot = p_ppc_mean, filename = "./img/ppc_mean.pdf")
