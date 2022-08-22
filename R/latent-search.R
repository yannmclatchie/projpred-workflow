library(glmnet)
library(cowplot)
library(ggplot2)
library(parallel)
library(LaplacesDemon)

# Poisson

data(PoissonExample)
x <- PoissonExample$x
y <- PoissonExample$y

fit_vanilla <- glmnet(x, y, family = "poisson")
fit_latent <- glmnet(x, log(y + 1), family = "gaussian")

plot(fit_vanilla, label = TRUE)
plot(fit_latent, label = TRUE)

cvfit_vanilla <- cv.glmnet(x, y, family = "poisson")
cvfit_latent <- cv.glmnet(x, log(y + 1), family = "gaussian")

plot(cvfit_vanilla)
plot(cvfit_latent)

