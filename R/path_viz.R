library(ggplot2)
library(dplyr)
library(latex2exp)

setwd("~/Desktop/projpred-workflow")

# read a csv
root <- paste0("./data/paths/thinning_paths.csv")
data <- read.csv(root)
data

p <- data %>%
  mutate(name=factor(name), rank=as.numeric(rank)) %>%
  group_by(name) %>%
  summarise(earliest_occurence = min(rank)) %>%
  right_join(data, by = "name") %>%
  filter(earliest_occurence < 25) %>%
  mutate(name1=fct_reorder2(name, rank, freq, .fun = weighted.mean, .desc=FALSE)) %>%
  ggplot(aes(x = name1, y = rank, fill = freq)) +
  geom_tile(aes(fill = freq)) + 
  scale_fill_gradient(low = "white", high = "steelblue", na.value = NA) +
  facet_grid(rho ~ thinning, scales = "free", labeller = "label_both") +
  coord_cartesian(expand = FALSE, ylim = c(0, 25)) +
  theme_classic() +
  labs(fill="Frequency") +
  ylab("Search path rank") +
  xlab("Covariate") +
  theme(
    strip.background = element_blank(), 
    strip.placement = "outside",
    axis.line=element_blank(),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_blank(),
    text=element_text(size=12), #change font size of all text
    axis.text=element_text(size=12), #change font size of axis text
    axis.text.x = element_blank(),#element_text(angle = 90),
    axis.title=element_text(size=12), #change font size of axis titles
    plot.title=element_text(size=12), #change font size of plot title
    legend.text=element_text(size=12), #change font size of legend text
    legend.title=element_text(size=12) #change font size of legend title   
  ) 
p

# save plot
file <- paste0("./img/thinning_path.pdf")
ggsave(plot = p, filename = file, dpi = "print")


