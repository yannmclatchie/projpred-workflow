# check if scam if installed
if (!require("scam")) install.packages("scam"); library("scam")

# load necessary libraries
library(dplyr)
library(brms)
library(projpred)
library(bayesplot)
library(parallel)
library(gamm4)
library(matrixStats)
library(geomtextpath)

# load the data
setwd("~/Desktop/projpred-workflow")
root <- "./data/benchmarks/student/"
d1=read.table(paste0(root, "student-mat.csv"),sep=";",header=TRUE)
d2=read.table(paste0(root, "student-por.csv"),sep=";",header=TRUE)
predictors <- c("school","sex","age","address","famsize","Pstatus","Medu","Fedu","traveltime","studytime","failures","schoolsup","famsup","paid","activities", "nursery", "higher", "internet", "romantic","famrel","freetime","goout","Dalc","Walc","health","absences","Mjob","Fjob","reason","guardian")
p <- length(predictors)
old_grades <- c("G1.x", "G2.x", "G3.x", "G1.y", "G2.y", "G3.y")
grades <- c("G1mat","G2mat","G3mat","G1por","G2por","G3por")
data <- merge(d1, d2, by=predictors) %>%
  rename_with(~ grades, all_of(old_grades))
data

# data pre-processing
data <- data %>%
  mutate(across(matches("G[1-3]..."), ~na_if(.,0))) %>%
  mutate(Gmat = rowMedians(as.matrix(dplyr::select(.,matches("G.mat"))), na.rm=TRUE),
         Gpor = rowMedians(as.matrix(dplyr::select(.,matches("G.por"))), na.rm=TRUE))
data_Gmat <- subset(data, is.finite(Gmat), select=c("Gmat",predictors))
data_Gpor <- subset(data, is.finite(Gpor), select=c("Gpor",predictors))

# remove those columns with only one unique value
data_Gmat <- data_Gmat %>% dplyr::select(where(~ n_distinct(.) > 1))
data_Gpor <- data_Gpor %>% dplyr::select(where(~ n_distinct(.) > 1))

(nmat <- nrow(data_Gmat))
(npor <- nrow(data_Gpor))

(p_Gmat <- dim(data_Gmat)[2] - 1)
(p_Gpor <- dim(data_Gpor)[2] - 1)

# fit the reference model
fit <- brm(
  formula = "Gmat ~ .", 
  family = "gaussian", 
  prior = set_prior(R2D2(mean_R2 = 0.25, prec_R2 = 50, cons_D2 = 1, autoscale = TRUE)), 
  data = data_Gmat
)

# diagnose the model computation
bayesplot::color_scheme_set("brightblue")
ratios <- neff_ratio(fit)
bayesplot::mcmc_neff(ratios)
rhats <- bayesplot::rhat(fit)
bayesplot::mcmc_rhat(rhats)

# diagnose the posterior predictive distribution
yrep <- brms::posterior_predict(fit, draws = 500)
bayesplot::ppc_dens_overlay(data_Gpor$Gpor, yrep[1:50, ])
bayesplot::ppc_hist(data_Gpor$Gpor, yrep[1:5, ])
bayesplot::ppc_stat(data_Gpor$Gpor, yrep, stat = "mean")
bayesplot::ppc_ecdf_overlay(data_Gpor$Gpor, yrep[1:50, ])

# perform projections
vs <- varsel(
  fit, 
  method = "forward", 
  nterms_max = p_Gpor
)

# compute submodel summaries
sel_df <- summary(vs, stats = c("elpd"))$selection

# fit different splines to the projpred results
spline <- gam(diff/diff.se ~ s(size, k = 10, m = 1), data = sel_df)
mono.spline <- scam(diff/diff.se ~ s(size, k = 10, bs = "mpi", m = 1), data = sel_df)

# update dataframe and plot
p <- sel_df %>%
  mutate(
    spline.diff = spline$fit*sel_df[,'diff.se'],
    spline.diff.se = sqrt(spline$sig2)*sel_df[,'diff.se'],
    mono.spline.diff = mono.spline$fit*sel_df[,'diff.se'],
    mono.spline.diff.se = sqrt(mono.spline$sig2)*sel_df[,'diff.se']
  ) %>%
  ggplot() +
    geom_pointrange(
      aes(
        x = size, 
        y = diff,
        ymin = diff - diff.se, 
        ymax = diff + diff.se
      )
    ) +
    geom_ribbon(
      aes(
        x = size,
        ymin = spline.diff - spline.diff.se, 
        ymax = spline.diff + spline.diff.se
      ), 
      alpha = 0.3, fill = "red"
    ) +
    geom_labelline(
      aes(x = size, y = spline.diff), 
      colour = "red", 
      label = "Spline", 
      hjust = "ymin"
    ) +
    geom_ribbon(
      aes(
        x = size,
        ymin = mono.spline.diff - mono.spline.diff.se, 
        ymax = mono.spline.diff + mono.spline.diff.se
      ), 
      alpha = 0.3, linetype = 2, fill = "blue"
    ) +
    geom_labelline(
      aes(x = size, y = mono.spline.diff), 
      colour = "blue", 
      label = "Monotonic spline", 
      hjust = "ymax"
    ) +
    geom_labelhline(
      aes(yintercept = 0),
      linetype = 2,
      colour = "grey",
      label = "Reference model",
      hjust = 0.1
    ) +
    ylab("ELPD difference") +
    xlab("Model Size") +
    theme_bw()
p
# save plot
#file <- paste0("./img/Gpor-spline.pdf")
file <- paste0("./img/Gmat-spline.pdf")
ggsave(plot = p, filename = file)

