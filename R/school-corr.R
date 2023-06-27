# load necessary libraries
if (!require("scam")) install.packages("scam"); library("scam")

library(dplyr)
library(purrr)
library(brms)
library(projpred)
library(matrixStats)
library(geomtextpath)
library(ggdendro)
library(patchwork)

# load the data
setwd("~/Desktop/projpred-workflow")
root <- "./data/benchmarks/student/"
d1=read.table(paste0(root, "student-mat.csv"),sep=";",header=TRUE)
d2=read.table(paste0(root, "student-por.csv"),sep=";",header=TRUE)
predictors <- c(
  "school",
  "sex",
  "age",
  "address",
  "famsize",
  "Pstatus",
  "Medu",
  "Fedu",
  "traveltime",
  "studytime",
  "failures",
  "schoolsup",
  "famsup",
  "paid",
  "activities", 
  "nursery", 
  #"higher", 
  "internet", 
  "romantic",
  "famrel",
  "freetime",
  "goout",
  "Dalc",
  "Walc",
  "health",
  "absences",
  "Mjob",
  "Fjob",
  "reason",
  "guardian"
)

p <- length(predictors)
old_grades <- c("G1.x", "G2.x", "G3.x", "G1.y", "G2.y", "G3.y")
grades <- c("G1mat","G2mat","G3mat","G1por","G2por","G3por")
data <- merge(d1, d2, by=predictors) %>%
  rename_with(~ grades, all_of(old_grades))

# data pre-processing
data <- data %>%
  mutate(across(matches("G[1-3]..."), ~na_if(.,0))) %>%
  mutate(
    Gmat = rowMedians(as.matrix(dplyr::select(.,matches("G.mat"))), na.rm=TRUE),
    Gpor = rowMedians(as.matrix(dplyr::select(.,matches("G.por"))), na.rm=TRUE))
data_Gmat <- subset(data, is.finite(Gmat), select=c("Gmat",predictors))
data_Gpor <- subset(data, is.finite(Gpor), select=c("Gpor",predictors))

# remove those columns with only one unique value
data_Gmat <- data_Gmat %>% dplyr::select(where(~ n_distinct(.) > 1))
data_Gpor <- data_Gpor %>% dplyr::select(where(~ n_distinct(.) > 1))

# convert character columns to factors
data_Gpor <- as.data.frame(unclass(data_Gpor), stringsAsFactors = TRUE)
data_Gmat <- as.data.frame(unclass(data_Gmat), stringsAsFactors = TRUE)

# standardise the data
datastd_Gmat <- data_Gmat
Gmatbin<-apply(data_Gmat[,predictors], 2, function(x) {length(unique(x))<10})
datastd_Gmat[,predictors[!Gmatbin]] <-scale(data_Gmat[,predictors[!Gmatbin]])
datastd_Gpor <- data_Gpor
Gporbin<-apply(data_Gpor[,predictors], 2, function(x) {length(unique(x))<10})
datastd_Gpor[,predictors[!Gporbin]] <-scale(data_Gpor[,predictors[!Gporbin]])

(nmat <- nrow(data_Gmat))
(npor <- nrow(data_Gpor))

(p_Gmat <- dim(data_Gmat)[2] - 1)
(p_Gpor <- dim(data_Gpor)[2] - 1)

# fit the reference model
fit <- brm(
  formula = Gmat ~ ., 
  family = gaussian, 
  prior = set_prior(R2D2(mean_R2 = 0.5, prec_R2 = 2, cons_D2 = 1)), 
  data = data_Gmat
)

################
## Projpred

# perform projections
vs <- cv_varsel(
  fit, 
  method = "forward", 
  nterms_max = p_Gmat,
  validate_search = TRUE
)

# compute submodel summaries
sel_df <- summary(vs, stats = c("elpd"))$selection

# KL divergence elbow
kl_df <- cbind(sel_df, vs$kl)
colnames(kl_df) <- c(colnames(sel_df), "kl")
( p_kl <- ggplot(kl_df, aes(x = size, y = kl)) +
  geom_point() +
  geom_line() + 
  ylab("KL divergence") +
  xlab("Model size") +
  theme_bw() +
  theme(legend.position="none") )
# save plot to tikz
#source("./R/aux/aux_plotting.R")
#save_tikz_plot(
#  p_kl, 
#  width = 3.5,
#  filename = "./tex/schools-kl.tex"
#)

# smoothed elpd differences
sel_df_not_null_mod <- filter(sel_df, size != 0)
mono.spline <- scam(
  diff/diff.se ~ s(size, k = 10, bs = "mpi", m = 2), 
  data = sel_df_not_null_mod
)
spline_df <- tibble(
  size = sel_df_not_null_mod$size[sel_df_not_null_mod$size>0],
  mono.spline.diff = mono.spline$fit*sel_df_not_null_mod[,'diff.se'],
  mono.spline.diff.se = sqrt(mono.spline$sig2)*sel_df_not_null_mod[,'diff.se']
)

################
## Dendrogram

dfun <- function(x) 1 - cor(x)
proj_and_pred <- function(x) {
  prj <- project(vs, solution_terms = x)
  preds <- proj_linpred(prj, transform = T, integrated = T)
  return (preds$pred)
}

cov_names <- names(data_Gmat)[names(data_Gmat) != "Gmat"]
single_cov_pred_list <- cov_names %>%
  map(proj_and_pred)
single_cov_pred_mat <- do.call("rbind", single_cov_pred_list) %>% t()
single_cov_pred_dists <- single_cov_pred_mat %>% dfun() %>% as.dist()
hc <- hclust(single_cov_pred_dists)
hc$labels <- cov_names
hcdata <- dendro_data(hc)

################
## Plotting
p_dendro <- ggplot() +
  geom_segment(
    data = segment(hcdata),
    aes(x = x, y = y, xend = xend, yend = yend)
  ) +
  scale_x_continuous(
    position = "top", breaks = hcdata$labels$x, labels = hcdata$labels$label
  ) +
  scale_y_reverse() +
  ylab("$dist{theta}{theta^prime}$") + 
  xlab("Covariate ($theta$)") + 
  coord_flip() +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()
        )

p_elpd_diff <- sel_df %>%
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
    data = spline_df,
    aes(
      x = size,
      ymin = mono.spline.diff - mono.spline.diff.se,
      ymax = mono.spline.diff + mono.spline.diff.se
    ),
    alpha = 0.2, linetype = 2, fill = "blue"
  ) +
  ylab("$delta elpd$") +
  xlab("Model size") +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
p_elpd_diff

p <- p_elpd_diff + p_dendro
p

# save plot to tikz
source("./R/aux/aux_plotting.R")
save_tikz_plot(p, width = 7, filename = "./tex/schools-dendro.tex")
