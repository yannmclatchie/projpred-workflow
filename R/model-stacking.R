# clean up environment
rm(list = ls())
gc(reset = TRUE)

# install libraries
if (!require("devtools")) {
  install.packages("devtools")
}

# install latest projpred dev branch from git
devtools::install_github("stan-dev/projpred", build_vignettes = TRUE)
library(projpred)
library(brms)
library(bayesplot)
library("Matrix")
library("MASS")

# set number of cores
options(mc.cores = parallel::detectCores())
set.seed(300416)

n <- 500
p <- 50

xi <- 0.34 # 0.59, 0.34, 0.25
rho <- 0.5 # 0, 0.5, 0.9

R1 <- matrix(rep(rho, 5 * 5), nrow=5)
R2 <- matrix(rep(rho, 5 * 5), nrow=5)
R3 <- matrix(rep(rho, 5 * 5), nrow=5)
R4 <- matrix(rep(rho, 35 * 35), nrow=35)
R <- matrix(bdiag(R1, R2, R3, R4), nrow=p)
diag(R) <- 1

w1 <- rep(xi, 5)
w2 <- rep(xi * 0.5, 5)
w3 <- rep(xi * 0.25, 5)
w4 <- rep(0, 35)
w <- c(w1, w2, w3, w4)

mu <- rep(0, p)

X <- mvrnorm(n=n, mu=mu, Sigma=R)
y <- X %*% w
data <- data.frame(y = y, X = X)

# fit the BRMS model
fit <- brm(formula="y ~ .", family="gaussian", prior = set_prior("horseshoe(3)"), data=data)

# using recommended parameters
ref_model <- cv_varsel(
  fit, 
  ndraws_pred = 100,
  validate_search = FALSE
)
summary <- ref_model$summary
summary["n"] <- n
summary["xi"] <- xi
summary["rho"] <- rho
summary["kl"] <- ref_model$kl
summary["loo_wts"] <- exp(summary["elpd.loo"]) / sum(exp(summary["elpd.loo"]))
dim(summary)
