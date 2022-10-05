library(ggplot2)
library(dplyr)
library(latex2exp)

setwd("~/Desktop/projpred-workflow")

## REGRESSION
## ----------

# concatenate the csvs
root <- "./data/pred_draws/regression.csv"
data <- read.csv(file = root)
data

# plot varying n and rho for different criteria
# ---------------------------------------------
data["ref.elpd"] <- data["elpd.loo"] - data["elpd.diff"]
data["ref.rmse"] <- data["rmse.loo"] - data["rmse.diff"]

# remove n = 50 case
data <- data %>%
  mutate_at(vars(ndraws_pred), factor)

# data wrangling
criteria_data <- data %>% 
  dplyr::select(c("n", "rho", "size", "rmse.loo", "ndraws_pred"))
criteria_data

hline_dat <- data %>% 
  group_by(n, rho) %>%
  summarise(
    rmse = mean(ref.rmse)
  )
hline_dat

size_4_data <- data %>%
  group_by(n, rho, ndraws_pred) %>%
  filter( (elpd.diff >= -4) ) %>%
  arrange(desc(as.numeric(size))) %>%
  slice(n = n()) %>%
  dplyr::select(size, rho, ndraws_pred, size)
size_4_data

# plot the results
p_ndraws_pred <- ggplot(data, aes(x = size, y = rmse.loo, colour = ndraws_pred)) +
  geom_ribbon(
    aes(
      ymin = rmse.loo - rmse.se, 
      ymax = rmse.loo + rmse.se, 
      fill = ndraws_pred, 
      colour = ndraws_pred
    ), 
    alpha = 0.2, linetype = 2
  ) +
  geom_line(size = 0.5) +
  geom_hline(
    data=hline_dat, 
    aes(yintercept=rmse), 
    linetype=2, 
    size=0.5, 
    colour="#000000"
  ) +
  geom_vline(
    data=size_4_data, 
    aes(xintercept=size, color=ndraws_pred, linetype=ndraws_pred), 
    size = 0.5
  ) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73")) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73")) +
  scale_linetype_manual(values=c("solid", "longdash", "dotted")) +
  coord_cartesian(clip = 'off') +
  ylab("RMSE") +
  xlab("Model size") +
  facet_grid(n ~ rho, scales = "free", labeller = "label_both") +
  theme_classic() +
  theme(
    strip.background = element_blank(), 
    strip.placement = "outside", 
    strip.text.x=element_text(size=16), 
    strip.text.y=element_text(size=16)
  ) +
  #theme(legend.position="none") + 
  theme(
    text=element_text(size=16), #change font size of all text
    axis.text=element_text(size=16), #change font size of axis text
    axis.title=element_text(size=16), #change font size of axis titles
    plot.title=element_text(size=16), #change font size of plot title
    legend.text=element_text(size=16), #change font size of legend text
    legend.title=element_text(size=16) #change font size of legend title   
  )
p_ndraws_pred

# save plot
file <- paste0("./img/ndraws_pred_rmse.pdf")
ggsave(plot = p_ndraws_pred, filename = file)


## CLASSIFICATION
## --------------

# concatenate the csvs
root <- "./data/pred_draws/classification.csv"
data <- read.csv(file = root)
data

# plot varying n and rho for different criteria
# ---------------------------------------------
data["ref.elpd"] <- data["elpd.loo"] - data["elpd.diff"]
data["ref.auc"] <- data["auc.loo"] - data["auc.diff"]

# remove n = 50 case
data <- data %>%
  mutate_at(vars(ndraws_pred), factor)

# data wrangling
criteria_data <- data %>% 
  dplyr::select(c("n", "rho", "size", "auc.loo", "ndraws_pred"))
criteria_data

hline_dat <- data %>% 
  group_by(n, rho) %>%
  summarise(
    auc = mean(ref.auc)
  )
hline_dat

size_4_data <- data %>%
  group_by(n, rho, ndraws_pred) %>%
  filter( (elpd.diff >= -4) ) %>%
  arrange(desc(as.numeric(size))) %>%
  slice(n = n()) %>%
  dplyr::select(size, rho, ndraws_pred, size)
size_4_data

# plot the results
p_ndraws_pred <- ggplot(data, aes(x = size, y = auc.loo, colour = ndraws_pred)) +
  geom_ribbon(
    aes(
      ymin = auc.loo - auc.se, 
      ymax = auc.loo + auc.se, 
      fill = ndraws_pred, 
      colour = ndraws_pred
    ), 
    alpha = 0.2, linetype = 2
  ) +
  geom_line(size = 0.5) +
  geom_hline(
    data=hline_dat, 
    aes(yintercept=auc), 
    linetype=2, 
    size=0.5, 
    colour="#000000"
  ) +
  geom_vline(
    data=size_4_data, 
    aes(xintercept=size, color=ndraws_pred, linetype=ndraws_pred), 
    size = 0.5
  ) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73")) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73")) +
  scale_linetype_manual(values=c("solid", "longdash", "dotted")) +
  coord_cartesian(clip = 'off') +
  ylab("AUC") +
  xlab("Model size") +
  facet_grid(n ~ rho, scales = "free", labeller = "label_both") +
  theme_classic() +
  theme(
    strip.background = element_blank(), 
    strip.placement = "outside", 
    strip.text.x=element_text(size=16), 
    strip.text.y=element_text(size=16)
  ) +
  #theme(legend.position="none") + 
  theme(
    text=element_text(size=16), #change font size of all text
    axis.text=element_text(size=16), #change font size of axis text
    axis.title=element_text(size=16), #change font size of axis titles
    plot.title=element_text(size=16), #change font size of plot title
    legend.text=element_text(size=16), #change font size of legend text
    legend.title=element_text(size=16) #change font size of legend title   
  )
p_ndraws_pred

# save plot
file <- paste0("./img/ndraws_pred_auc.pdf")
ggsave(plot = p_ndraws_pred, filename = file)

