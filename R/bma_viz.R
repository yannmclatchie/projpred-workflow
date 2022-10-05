library(ggplot2)
library(dplyr)
library(latex2exp)

setwd("~/Desktop/projpred-workflow")

## BMA-proj
## ----------

# concatenate the csvs
proj_root <- "./data/bma/bma_proj.csv"
proj_data <- read.csv(file = proj_root)
proj_data

# plot varying n and rho for different criteria
# ---------------------------------------------

# plot the criteria
p_proj <- ggplot(proj_data, aes(x = size, y = mlpd.diff, colour = data)) +
  geom_line() + 
  geom_vline(xintercept = 15, linetype=2, size=0.35, color='grey') +
  facet_grid(n ~ rho, scales = "free", labeller = label_both) +
  ylab("MLPD difference") +
  xlab("Model size") +
  theme_classic()
p_proj

# save the plots
ggsave(plot = p_proj, filename = "./img/bma_proj.pdf")

## BMA-ref
## ----------

# concatenate the csvs
ref_root <- "./data/bma/bma_ref.csv"
ref_data <- read.csv(file = ref_root)
ref_data

# preprocessing
ref_data["mlpd"] <- ref_data["elpd"] / ref_data["n"]

ref_mlpd_data <- ref_data %>% 
  group_by(data, n, rho) %>% 
  filter(size == max(size)) %>%
  mutate(ref_mlpd = mlpd) %>%
  dplyr::select(data, n, rho, ref_mlpd) 
plot_data <- left_join(ref_data, ref_mlpd_data, by = c("data", "n", "rho"))
plot_data["mlpd.diff"] <- plot_data["mlpd"] - plot_data["ref_mlpd"]

# plot varying n and rho for different criteria
# ---------------------------------------------

# plot the criteria
p_ref <- ggplot(plot_data, aes(x = size, y = mlpd.diff, colour = data)) +
  geom_line() + 
  geom_vline(xintercept = 15, linetype=2, size=0.35, color='grey') +
  facet_grid(n ~ rho, scales = "free", labeller = label_both) +
  ylab("MLPD difference") +
  xlab("Model size") +
  theme_classic()
p_ref

# save the plots
ggsave(plot = p_ref, filename = "./img/bma_ref.pdf")
