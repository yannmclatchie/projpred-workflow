library(ggplot2)
library(projpred)
library(dplyr)
library(tidyr)
library(ggdendro)
source('R/aux/aux_plotting.R')
### ELPD plot
load('R/dat/case_study_4.1_500.RData')
elpd_summary <- bind_rows(list(forward_loo=vs_forward_loo$summary,
                               L1_loo=vs_L1_loo$summary,
                               forward_loo_validated=vs_forward_loo_validated$summary,
                               forward_kfold=rename(vs_forward_kfold$summary,elpd.loo=elpd.kfold)),.id='varsel') %>%
  mutate(varsel_lab=case_when(varsel=='forward_loo' ~ 'Single-path forward selection',
                              varsel=='L1_loo' ~ 'LOO L1 selection',
                              varsel=='forward_loo_validated' ~ 'Cross-validated forward selection',
                              TRUE ~ 'K-fold forward selection'))
max_size_displayed <- 30
elpd_ref <- filter(elpd_summary,varsel=='forward_loo_validated') %>% slice(1) %>% mutate(elpd_ref=elpd.loo-diff) %>% .$elpd_ref
mod_size_candidates <- c('confidence_heuristic' = with(vs_forward_loo_validated$summary,size[min(which(diff + diff.se>=0))]),
                         'elpd_difference_heuristic'=spline_df$size[min(which(spline_df$mono.spline.diff >= -4))])
elpd_plot <- filter(elpd_summary,size<=max_size_displayed,varsel %in% c('forward_loo','forward_loo_validated')) %>%
  ggplot(aes(x=size,
             y=elpd.loo,
             ymin=elpd.loo-se,
             ymax=elpd.loo+se,
             col=varsel_lab)) +
  geom_hline(yintercept=elpd_ref,colour = "red", linetype = "longdash") +
  geom_vline(xintercept = 9, colour = "grey", linetype = "longdash") +
  geom_vline(xintercept = 6, colour = "grey", linetype = "longdash") +
  geom_pointrange(position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0)) +
  geom_line(position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0)) +
  annotate("text", x = 1.25, y = elpd_ref + 10, colour = "red", label = "Reference model elpd",size=3) +
  scale_color_manual(values=c('black','grey'),name='') +
  scale_x_continuous(breaks = 0:max_size_displayed) +
  xlab('Model size') +
  ylab('elpd') +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = c(0.8,0.55))
save_tikz_plot(elpd_plot,width=6.5,height=4,filename = 'tex/elpd_plot_case_study_4.1.tex')



suggest_size(vs_forward_loo,thres_elpd = -4)
suggest_size(vs_forward_loo_validated,thres_elpd = -4)

suggest_size(vs_forward_loo)
suggest_size(vs_forward_loo_validated)
### DENDROGRAM
dfun <- function(x) 1 - cor(x)
# code from https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/
# to enable coloring specific segments of dendrogram
dendro_data_k <- function(hc, k) {

  hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
  seg       <-  hcdata$segments
  labclust  <-  cutree(hc, k)[hc$order]
  segclust  <-  rep(0L, nrow(seg))
  heights   <-  sort(hc$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)

  for (i in 1:k) {
    xi      <-  hcdata$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }

  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  hcdata$segments$clust  <-  segclust
  hcdata$segments$line   <-  as.integer(segclust < 1L)
  hcdata$labels$clust    <-  labclust

  hcdata
}

cov_names <- unique(c(paste0('X.',1:100),vs_forward_loo_validated$solution_terms))
single_cov_preds <- lapply(cov_names,function(x){
  prj <- project(mod_ref, solution_terms = x)
  preds <- proj_linpred(prj, transform = T, integrated = T)
  preds$pred
}) %>% do.call('rbind',.) %>%
  t()
single_cov_preds_dist <- dfun(single_cov_preds) %>% as.dist()

hc <- hclust(single_cov_preds_dist)
hc$labels <- cov_names
hcdata <- dendro_data_k(hc,k=30)
relevant_clusters <- unique(hcdata$labels$clust[hcdata$labels$label %in% paste0('X.',1:15)])
hcdata$segments$clust <- hcdata$segments$clust %in% relevant_clusters
hcdata$labels$clust <- hcdata$labels$clust %in% relevant_clusters
cov_selected <- vs_forward_loo_validated$solution_terms[seq_len(6)]
hcdata$labels$label <- ifelse(hcdata$labels$label %in% c('X.9'),paste0('\n',hcdata$labels$label),hcdata$labels$label)
cov_selected <- ifelse(cov_selected%in% c('X.9'),paste0('\n',cov_selected),cov_selected)
p_dendro <- ggplot() +
  geom_segment(
    data = segment(hcdata),
    aes(x = x, y = y, xend = xend, yend = yend,col=clust)
  ) +
  scale_x_continuous(
    position = "top", breaks = hcdata$labels$x[hcdata$labels$label %in% cov_selected], labels = hcdata$labels$label[hcdata$labels$label %in% cov_selected]
  ) +
  scale_y_reverse() +
  scale_color_manual(values=c('black','blue')) +
  ylab("$dist{theta}{theta^prime}$") +
  xlab("Covariate ($theta$)") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = 'none'
  )
save_tikz_plot(p_dendro,width=3,height=4,filename = 'tex/dendrogram_case_study_4.1.tex')
