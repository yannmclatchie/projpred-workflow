# load necessary libraries
library(tidyverse)
library(brms)
library(projpred)
library(parallel)

# load data generating scripts
source("./R/generate.R")

# define experiment
experiment <- function(rho, regul){
  n <- 200
  
  # simulate data
  data <- .generate(rho = rho, n = n)
  
  # fit the BRMS model
  fit <- brm(
    formula = "y ~ .", 
    family = "gaussian", 
    prior = set_prior("horseshoe(11)"), 
    data = data
  )
  
  # perform projection predictive inference
  vs <- varsel(
    fit,
    method = "L1", 
    regul = regul,
    nterms_max = 100
  )
  paths <- vs$solution_terms
  
  for (i in 1:1) {
    # simulate data
    data <- .generate(rho, n)
    
    # fit the BRMS model
    fit <- brm(
      formula = "y ~ .", 
      family = "gaussian", 
      prior = set_prior("horseshoe(11)"), 
      data = data
    )
    
    # perform projection predictive inference
    vs <- varsel(
      fit,
      method = "L1", 
      regul = regul,
      nterms_max = 100
    )
    paths <- cbind(paths, vs$solution_terms)
  }
  
  # helper function to compute occurrences of covariate rankings
  .count_occs <- function(x) {
    occs <- which(paths == x, arr.ind = TRUE)
    while (nrow(occs) < ncol(paths)) {
      pad <- data.frame(row = 0, col = 0)
      occs <- rbind(occs, pad)
    }
    occs[occs == 0] <- NA
    res <- as.data.frame(table(occs[,"row"]))
    res$Freq <- res$Freq / sum(res$Freq)
    return (res)
  }
  
  els <- unique(as.vector(as.matrix(paths)))
  occs <- purrr::map(els, .count_occs)
  names(occs) <- els
  res <- do.call("rbind", occs)
  
  res_names <- sub("(.*?..*?).(.*?)", "\\1", rownames(res))
  res_names <- sub("\\..*", "", res_names)
  res["name"] <- res_names
  res["rho"] <- rho
  res["regul"] <- regul
  res <- rename(res, rank=Var1, freq=Freq)
  return(res)
}

# vary both the number of data observations and the correlation between parameters
rhos <- c(0.5)
reguls <- c(0.5, 1 - 1e-04)
options <- expand.grid(rho=rhos, regul=reguls)

# perform experiment
res <- parallel::mcMap(
  f = experiment, 
  rho = options[, "rho"], 
  regul = options[, "regul"], 
  mc.cores = 1
)
# concatenate experiment results
data <- do.call("rbind", res)

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
  facet_grid(rho ~ regul, scales = "free", labeller = "label_both") +
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


# write the table
setwd("/scratch/work/mclatcy1/projpred-workflow")
csv_name <- paste0("./data/paths/regul/regul_paths.csv")
ff <- file(csv_name, open="w")
write.table(df, file = ff, sep = ",", row.names = FALSE)
close(ff)
