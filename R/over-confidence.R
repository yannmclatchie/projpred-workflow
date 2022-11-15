library(simstudy)
library(bayesplot)

# define study
def <- defData(varname = "x1", formula = "0", variance = 1,
               dist = "normal")
def <- defData(def, varname = "x2", formula = "0", variance = 1,
               dist = "normal")
def <- defRepeat(nVars = 97, prefix = "x", formula = "0", variance = 1,
                 dist = "normal")
def <- defData(def, varname = "x100", formula = "x2", variance = 0.1,
                dist = "normal")

def <- defData(def, "y", formula = "x1 + x2", variance = 1,
                dist = "normal")

# set seed
set.seed(87261)

# generate data from the study
dd <- genData(1000, def)

# fit a reference model to the data
fit <- brm(
  formula = "y ~ .", 
  family = "gaussian", 
  #prior = set_prior(R2D2(mean_R2 = 0.3, prec_R2 = 3)), 
  data = dd
)

# refit the oracle submodel
refit <- brm(
  formula = "y ~ x1 + x2", 
  family = "gaussian", 
  #prior = set_prior(R2D2(mean_R2 = 0.3, prec_R2 = 3)), 
  data = dd
)

# project the model onto the true DGP
vs <- varsel(fit, nterms_max = 10)
prj <- project(
  vs, 
  solution_terms = c("x1", "x2")
)
prj_mat <- as.matrix(prj)

# plot the three posteriors
p <- ggplot() + 
  mcmc_areas(prj_mat) + 
  mcmc_areas(fit, pars = c("b_Intercept", "b_x1", "b_x2", "sigma")) +
  mcmc_areas(refit, pars = c("b_Intercept", "b_x1", "b_x2", "sigma"))
p

# plot the marginals of the correlated covariates
mcmc_scatter(fit, pars = c("b_x2", "b_x100"),
             size = 1.5, alpha = 0.5)

