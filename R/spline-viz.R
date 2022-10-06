library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/Desktop/projpred-workflow")
root <- "./data/smoothed/spline_stability.csv"
df <- read.csv(root)

## Model size selection
## --------------------
loo.elpd.delta.size <- df %>% group_by(n, rho, run) %>%
  filter( (diff >= -4) ) %>%
  arrange(desc(as.numeric(size))) %>%
  slice(n = n()) %>%
  dplyr::select(size)
spline.elpd.delta.size <- df %>% group_by(n, rho, run) %>%
  filter( (spline.elpd.diff >= -4) ) %>%
  arrange(desc(as.numeric(size))) %>%
  slice(n = n()) %>%
  dplyr::select(size)
mono.spline.elpd.delta.size <- df %>% group_by(n, rho, run) %>%
  filter( (mono.spline.elpd.diff >= -4) ) %>%
  arrange(desc(as.numeric(size))) %>%
  slice(n = n()) %>%
  dplyr::select(size)
loo.elpd.ci.size <- df %>% group_by(n, rho, run) %>%
  filter( (elpd + se >= ref.loo.elpd) ) %>%
  arrange(desc(as.numeric(size))) %>%
  slice(n = n()) %>%
  dplyr::select(size)
spline.elpd.ci.size <- df %>% group_by(n, rho, run) %>%
  filter( (spline.elpd + spline.elpd.se >= ref.loo.elpd) ) %>%
  arrange(desc(as.numeric(size))) %>%
  slice(n = n()) %>%
  dplyr::select(size)
mono.spline.elpd.ci.size <- df %>% group_by(n, rho, run) %>%
  filter( (mono.spline.elpd + mono.spline.elpd.se >= ref.loo.elpd) ) %>%
  arrange(desc(as.numeric(size))) %>%
  slice(n = n()) %>%
  dplyr::select(size)
loo.elpd.switch.size <- min(loo.elpd.ci.size, loo.elpd.delta.size)
spline.elpd.switch.size <- min(spline.elpd.ci.size, spline.elpd.delta.size)
mono.spline.elpd.switch.size <- min(
  mono.spline.elpd.ci.size, mono.spline.elpd.delta.size
)

sel_df <- data.frame(
    n = loo.elpd.delta.size$n,
    rho = loo.elpd.delta.size$rho,
    run = loo.elpd.delta.size$run,
    loo.elpd.delta = loo.elpd.delta.size$size,
    spline.elpd.delta = spline.elpd.delta.size$size,
    mono.spline.elpd.delta = mono.spline.elpd.delta.size$size,
    loo.elpd.ci = loo.elpd.ci.size$size,
    spline.elpd.ci = spline.elpd.ci.size$size,
    mono.spline.elpd.ci = mono.spline.elpd.ci.size$size
  ) %>% rowwise() %>%
  mutate(
    loo.elpd.switch.size = min(loo.elpd.ci, loo.elpd.delta),
    spline.elpd.switch.size = min(spline.elpd.ci, spline.elpd.delta),
    mono.spline.elpd.switch.size = min(
      mono.spline.elpd.ci, mono.spline.elpd.delta
    )
  )

## Plotting
## --------

# model selection stability
p_sel_stability <- sel_df %>%
  reshape2::melt(id.vars = c("n", "rho", "run")) %>%
  group_by(n, rho, variable) %>% 
  mutate(
    mean_size = mean(value),
    lower_size = quantile(value, probs = 0.32),
    upper_size = quantile(value, probs = 0.68)
  ) %>%
  dplyr::select(n, rho, variable, mean_size, lower_size, upper_size) %>%
  distinct() %>%
  ggplot() +
  geom_pointrange(
    aes(
      x = variable, 
      y = mean_size,
      ymin = lower_size, 
      ymax = upper_size, 
      colour = variable
    )
  ) +
  labs(colour="Procedure") +
  facet_grid(n ~ rho, scales = "free", labeller = "label_both") +
  ylab("Size") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    axis.title.x=element_blank(), 
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) 
p_sel_stability
ggsave(
  plot = p_sel_stability, 
  filename = "./img/spline_selection_stability.pdf"
)

# ELPDs
means = df %>% 
  group_by(n, rho, size) %>% 
  summarize(
    mean_loo = mean(elpd),
    mean_spline = mean(spline.elpd),
    mean_mono_spline = mean(mono.spline.elpd)
  )

p_smoothed_elpds <- df %>%
  ggplot() + 
    geom_line(aes(x = size, y = elpd, group = run, color = 'LOO'), 
              size = 0.5, alpha = 0.15) +
    geom_line(aes(x = size, y = spline.elpd, group = run, color = 'Spline'), 
              size = 0.5, alpha = 0.15) +
    geom_line(aes(x = size, y = mono.spline.elpd, group = run, color = 'Mono'), 
              size = 0.5, alpha = 0.15) +
    geom_line(data = means, aes(x = size, y = mean_loo, color = 'LOO'), size = 1) +
    geom_line(data = means, aes(x = size, y = mean_spline, color = 'Spline'), size = 1) +
    geom_line(data = means, aes(x = size, y = mean_mono_spline, color = 'Mono'), size = 1) +
    facet_grid(n ~ rho, scales = "free", labeller = "label_both") +
    xlab("Size") +
    ylab("ELPD") +
    theme_bw() 
p_smoothed_elpds
ggsave(
  plot = p_smoothed_elpds, 
  filename = "./img/smoothed_elpd.pdf"
)

# ELPD diffs
means = df %>% 
  group_by(n, rho, size) %>% 
  summarize(
    mean_loo = mean(diff),
    mean_spline = mean(spline.elpd.diff),
    mean_mono_spline = mean(mono.spline.elpd.diff)
  )

p_smoothed_elpd_diffs <- df %>%
  ggplot() + 
  geom_line(aes(x = size, y = diff, group = run, color = 'LOO'), 
            size = 0.5, alpha = 0.15) +
  geom_line(aes(x = size, y = spline.elpd.diff, group = run, color = 'Spline'), 
            size = 0.5, alpha = 0.15) +
  geom_line(aes(x = size, y = mono.spline.elpd.diff, group = run, color = 'Mono'), 
            size = 0.5, alpha = 0.15) +
  geom_line(data = means, aes(x = size, y = mean_loo, color = 'LOO'), size = 1) +
  geom_line(data = means, aes(x = size, y = mean_spline, color = 'Spline'), size = 1) +
  geom_line(data = means, aes(x = size, y = mean_mono_spline, color = 'Mono'), size = 1) +
  facet_grid(n ~ rho, scales = "free", labeller = "label_both") +
  xlab("Size") +
  ylab("ELPD difference") +
  theme_bw() 
p_smoothed_elpd_diffs
ggsave(
  plot = p_smoothed_elpd_diffs, 
  filename = "./img/smoothed_elpd_diff.pdf"
)

# Selected model ELPDs
sel_df_long <- pivot_longer(
  sel_df,
  cols = - one_of(c('n','rho','run')),
  names_to='Procedure',
  values_to='size'
)

sel_elpd_df <- inner_join(sel_df_long, df, by=c('n','rho','run','size')) %>%
  select(n, rho, run, Procedure, elpd, se, diff, diff.se, ref.loo.elpd)

p_sel_elpd <- sel_elpd_df %>% group_by(n, rho, Procedure) %>%
  summarize(
    mean_elpd = mean(elpd),
    mean_se = mean(se)
  ) %>%
  ggplot() +
  geom_pointrange(
    aes(
      x = Procedure, 
      y = mean_elpd,
      ymin = mean_elpd - mean_se, 
      ymax = mean_elpd + mean_se, 
      colour = Procedure
    )
  ) +
  labs(colour="Procedure") +
  facet_grid(n ~ rho, scales = "free", labeller = "label_both") +
  ylab("ELPD") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    axis.title.x=element_blank(), 
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) 
p_sel_elpd
ggsave(
  plot = p_sel_elpd, 
  filename = "./img/sel_elpd.pdf"
)

p_sel_elpd_diff <- sel_elpd_df %>% group_by(n, rho, Procedure) %>%
  summarize(
    mean_diff = mean(diff),
    mean_diff_se = mean(diff.se)
  ) %>%
  ggplot() +
  geom_pointrange(
    aes(
      x = Procedure, 
      y = mean_diff,
      ymin = mean_diff - mean_diff_se, 
      ymax = mean_diff + mean_diff_se, 
      colour = Procedure
    )
  ) +
  labs(colour="Procedure") +
  facet_grid(n ~ rho, scales = "free", labeller = "label_both") +
  ylab("ELPD difference") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    axis.title.x=element_blank(), 
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) 
p_sel_elpd_diff
ggsave(
  plot = p_sel_elpd_diff, 
  filename = "./img/sel_elpd_diff.pdf"
)
