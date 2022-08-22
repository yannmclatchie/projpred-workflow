library(vroom)
library(ggplot2)
library(dplyr)
library(latex2exp)
library(geomtextpath)

# concatenate the csvs
root <- "./data/elbows/"
files <- paste0(root, list.files(path = root))
data <- vroom(files, delim = ",")

# plot varying n and rho for different criteria
# ---------------------------------------------
data["ref.elpd"] <- data["elpd"] - data["elpd.diff"]
data["ref.rmse"] <- data["rmse"] - data["rmse.diff"]
hline_dat <- data %>%
  group_by(n, rho) %>%
  summarise(
    ref.elpd = mean(ref.elpd),
    ref.rmse = mean(ref.rmse)
  )

# plot the criteria
p_loo <- ggplot(data, aes(x = size, y = loo_wts)) +
  geom_point() +
  geom_smooth() + 
  geom_vline(xintercept = 15, linetype=2, size=0.35, color='grey') +
  facet_grid(n ~ rho, scales = "free", labeller = label_both) +
  ylab("LOO weights") +
  xlab("Model size") +
  theme_classic() +
  theme(legend.position="none")
p_kl <- ggplot(data, aes(x = size, y = kl)) +
  geom_point() +
  geom_smooth() + 
  geom_vline(xintercept = 15, linetype=2, size=0.35, color='grey') +
  geom_hline(yintercept = 0, linetype=3, size=0.35, color='salmon') +
  facet_grid(n ~ rho, scales = "free", labeller = label_both) +
  ylab("KL divergence") +
  xlab("Model size") +
  theme_classic() +
  theme(legend.position="none")
p_rmse <- ggplot(data, aes(x = size, y = rmse)) +
  geom_point() +
  geom_smooth() + 
  geom_hline(
    data=hline_dat, 
    aes(yintercept=ref.rmse), 
    linetype=3, 
    size=0.35, 
    colour="salmon"
  ) +
  geom_vline(xintercept = 15, linetype=2, size=0.35, color='grey') +
  facet_grid(n ~ rho, scales = "free", labeller = label_both) +
  ylab("RMSE") +
  xlab("Model size") +
  theme_classic() +
  theme(legend.position="none")
p_elpd <- ggplot(data, aes(x = size, y = elpd)) +
  geom_point() +
  geom_smooth() + 
  geom_hline(
    data=hline_dat, 
    aes(yintercept=ref.elpd), 
    linetype=3, 
    size=0.35, 
    colour="salmon"
  ) +
  geom_vline(xintercept = 15, linetype=2, size=0.35, color='grey') +
  facet_grid(n ~ rho, scales = "free", labeller = label_both) +
  ylab("ELPD") +
  xlab("Model size") +
  theme_classic() +
  theme(legend.position="none")

# save the plots
ggsave(plot = p_loo, filename = "./img/loo_wts.pdf", width=5, height=5)
ggsave(plot = p_kl, filename = "./img/kl.pdf", width=5, height=5)
ggsave(plot = p_elpd, filename = "./img/elpd.pdf", width=5, height=5)
ggsave(plot = p_rmse, filename = "./img/rmse.pdf", width=5, height=5)

# plot fixed n criteria
# ---------------------

Ns <- c(50, 100, 200, 500)

for (N in Ns){
  print(paste0("Running visualisations for ",N," observations ..."))
  
  # data wrangling
  fixed_n_data <- data %>% 
    filter(n == N) %>% 
    select(c("size", "rho", "elpd.diff", "kl",  "rmse"))
  
  plot_data <- reshape2::melt(fixed_n_data, id.vars = c("size", "rho"))
  
  hline_dat <- data %>% 
    filter(n == N) %>% 
    summarise(
      elpd.diff = 4,
      rmse = mean(ref.rmse),
      kl = 0
    ) %>%
    reshape2::melt()
  names(hline_dat) <- c("criterion", "base")
  
  plot_data %>%
    group_by(rho, variable) %>%
    filter( (variable == "elpd.diff") ) %>%
    filter( (value <= 4) ) %>%
    arrange(desc(as.numeric(size))) %>%
    slice(n = n())
  names(plot_data) <- c("size", "rho", "criterion", "value")
  
  size_3_data <- plot_data %>%
    group_by(rho, criterion) %>%
    filter( (criterion == "elpd.diff") ) %>%
    filter( (value <= 3) ) %>%
    arrange(desc(as.numeric(size))) %>%
    slice(n = n()) %>%
    select(criterion, rho, size) %>%
    slice(rep(1, length(unique(plot_data$criterion)))) %>%
    mutate(
      criterion = rep(unique(plot_data$criterion))
    )
  
  size_4_data <- plot_data %>%
    group_by(rho, criterion) %>%
    filter( (criterion == "elpd.diff") ) %>%
    filter( (value <= 4) ) %>%
    arrange(desc(as.numeric(size))) %>%
    slice(n = n()) %>%
    select(criterion, rho, size) %>%
    slice(rep(1, length(unique(plot_data$criterion)))) %>%
    mutate(
      criterion = rep(unique(plot_data$criterion))
    )
  
  size_5_data <- plot_data %>%
    group_by(rho, criterion) %>%
    filter( (criterion == "elpd.diff") ) %>%
    filter( (value <= 5) ) %>%
    arrange(desc(as.numeric(size))) %>%
    slice(n = n()) %>%
    select(criterion, rho, size) %>%
    slice(rep(1, length(unique(plot_data$criterion)))) %>%
    mutate(
      criterion = rep(unique(plot_data$criterion))
    )
  
  print("Plotting ...")
  # plot the results
  p_criteria <- ggplot(plot_data, aes(x = size, y = value)) +
    geom_point() +
    geom_smooth() + 
    geom_hline(
      data=hline_dat, 
      aes(yintercept=base), 
      linetype=3, 
      size=0.35, 
      colour="salmon"
    ) +
    geom_textvline(
      data=size_3_data, 
      aes(xintercept=size), 
      label = TeX("$\\Delta = 3$"),
      # size = 0.5, 
      linetype = 2,
      vjust = 1.5, 
      linewidth = 0.35, 
      colour="seagreen2"
    ) +
    geom_textvline(
      data=size_4_data, 
      aes(xintercept=size), 
      label = TeX("$\\Delta = 4$"),
      # size = 0.5,
      linetype=2, 
      vjust = 0, 
      linewidth = 0.35,
      colour="seagreen3"
    ) +
    geom_textvline(
      data=size_5_data, 
      aes(xintercept=size), 
      label = TeX("$\\Delta = 5$"),
      # size = 0.5,
      linetype = 2, 
      vjust = -1.5, 
      linewidth = 0.35, 
      colour="seagreen4"
    ) +
    coord_cartesian(clip = 'off') +
    ylab("Criterion") +
    xlab("Model size") +
    facet_grid(criterion ~ rho, scales = "free") +
    theme_classic() +
    theme(legend.position="none")
  
  # save the plot
  print("Saving ...")
  file <- paste0("./img/concrete-regression/",N,".pdf")
  ggsave(plot = p_criteria, filename = file, width=5, height=5)
  print("Done!")
}

