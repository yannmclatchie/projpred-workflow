library(vroom)
library(ggplot2)
library(dplyr)
library(latex2exp)

setwd("~/Desktop/projpred-workflow")

# concatenate the csvs
root <- "./data/method/"
files <- paste0(root, list.files(path = root))
data <- vroom(files, delim = ",")

# plot varying n and rho for different criteria
# ---------------------------------------------
data["ref.elpd"] <- data["elpd"] - data["elpd.diff"]
data["ref.rmse"] <- data["rmse"] - data["rmse.diff"]

# remove n = 50 case
data <- data %>%
  filter( (n > 50) )

# data wrangling
criteria_data <- data %>% 
  dplyr::select(c("n", "rho", "size", "rmse", "method"))

hline_dat <- data %>% 
  group_by(n, rho) %>%
  summarise(
    rmse = mean(ref.rmse)
  )
hline_dat

size_4_data <- data %>%
  group_by(n, rho, method) %>%
  filter( (elpd.diff >= -4) ) %>%
  arrange(desc(as.numeric(size))) %>%
  slice(n = n()) %>%
  dplyr::select(size, rho, method, size)
size_4_data

# plot the results
p_method <- ggplot(data, aes(x = size, y = rmse, colour = method)) +
  geom_ribbon(
    aes(ymin = rmse - rmse.se, ymax = rmse + rmse.se, fill = method, colour = method), 
    alpha = 0.3, linetype = 2
  ) +
  geom_hline(
    data=hline_dat, 
    aes(yintercept=rmse), 
    linetype=2, 
    size=0.5, 
    colour="#000000"
  ) +
  #geom_vline(
  #  data=size_4_data, 
  #  aes(xintercept=size, color=method, linetype=method), 
  #  size = 0.5
  #) +
  geom_line(size = 0.5) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  facet_grid(n ~ rho, scales = "free", labeller = "label_both") +
  coord_cartesian(clip = 'off') +
  ylab("RMSE") +
  xlab("Model size") +
  theme_classic() +
  theme(
    legend.position="none",
    strip.background = element_blank(), 
    strip.placement = "outside", 
    strip.text.x=element_text(size=16), 
    strip.text.y=element_text(size=16),
    text=element_text(size=16), #change font size of all text
    axis.text=element_text(size=16), #change font size of axis text
    axis.title=element_text(size=16), #change font size of axis titles
    plot.title=element_text(size=16), #change font size of plot title
    legend.text=element_text(size=16), #change font size of legend text
    legend.title=element_text(size=16) #change font size of legend title   
  )
p_method

# save plot
file <- paste0("./img/method.pdf")
ggsave(plot = p_method, filename = file)
