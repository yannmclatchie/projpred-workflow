#library(ggplot2)
#library(dplyr)
#library(tidyr)
library(brms)
library(projpred)
#library(ggdendro)
source('R/generate.R')
#source('R/aux/aux_plotting.R')
# move two folders up from sbatch folder
.libPaths('../R/x86_64-pc-linux-gnu-library/4.1/')

set.seed(4321)
rho <- 0.9
n <- 500
dat <- .generate(rho=rho,n=n)
prior_ref <- prior(R2D2(mean_R2=0.3,prec_R2=5,cons_D2=1))
mod_ref <- brm(as.formula(paste0('y~',paste(paste0('X.',1:100),collapse='+'))),data=dat,prior=prior_ref,family=gaussian(),refresh=0,iter=2000,cores=4)

# rhats <- rhat(mod_ref)
# loo_ref <- loo(mod_ref)
# neff_ratio(mod_ref)

vs_forward_loo <- cv_varsel(mod_ref,cv_method='LOO',method = "forward",validate_search=F,nterms_max = 30)
vs_L1_loo <- cv_varsel(mod_ref,cv_method='LOO',method = "L1",validate_search=F,nterms_max = 30)
vs_forward_kfold <- cv_varsel(mod_ref,cv_method='kfold',method = "forward",nterms_max = 30)
size_maximum_elpd <- vs_forward_loo$summary$size[which.max(vs_forward_loo$summary$elpd.loo)]
print(paste('running cross validation of search path with max terms =',size_maximum_elpd))
vs_forward_loo_validated <- cv_varsel(mod_ref,cv_method='LOO',method = "forward",validate_search=T,nterms_max = size_maximum_elpd)
save(vs_forward_loo,vs_L1_loo,vs_forward_kfold,file=paste0('R/dat/varsel_objects_',n,'.RData'))


# elpd_summary <- bind_rows(list(forward_loo=vs_forward_loo$summary,
#                           L1_loo=vs_L1_loo$summary,
#                           forward_kfold=rename(vs_forward_kfold$summary,elpd.loo=elpd.kfold)),.id='varsel')
#
# ggplot(filter(elpd_summary,size<=20)) +
# geom_point(aes(size,elpd.loo,col=varsel)) +
# geom_line(aes(size,elpd.loo,col=varsel)) +
# geom_errorbar(aes(x=size,ymin=elpd.loo-se,ymax=elpd.loo+se,col=varsel),width=0.2) +
# geom_hline(yintercept=elpd_summary$elpd.loo[1] - elpd_summary$diff[1],linetype='dashed') +
# theme_bw()
#
#
#
# #varsel_forward_loo <- varsel(mod_ref,method = "forward",nterms_max = 30)
#
# # ggsave(paste0('R/fig/cv_varsel_forward_loo_',,'.png',plot(vs_forward_loo))
# # ggsave('R/fig/cv_varsel_forward_kfold.png',plot(vs_forward_kfold))
# # ggsave('R/fig/cv_varsel_L1_loo.png',plot(vs_forward))
# # ggsave('R/fig/cv_varsel_forward_loo_validated.png',plot(vs_forward_loo_validated))
#
# ggplot(tibble(model_size=seq(0,length(vs_forward$kl)-1),KL=vs_forward$kl),aes(model_size,KL)) +
# geom_point() +
# geom_line() +
# xlab('Model size') +
# ylab('Divergence from the reference model') +
# theme_bw()
#
# suggest_size(vs_forward_kfold,thres_elpd = -4)
# suggest_size(vs_forward_kfold)
#
# dfun <- function(x) 1 - cor(x)
#
# cov_names <- unique(c(paste0('X.',1:100),vs_forward_loo_validated$solution_terms))
# single_cov_preds <- lapply(cov_names,function(x){
#                       prj <- project(mod_ref, solution_terms = x)
#                       preds <- proj_linpred(prj, transform = T, integrated = T)
#                       preds$pred
#                     }) %>% do.call('rbind',.) %>%
#                     t()
# single_cov_preds_dist <- dfun(single_cov_preds) %>% as.dist()
#
# hc <- hclust(single_cov_preds_dist)
# hc$labels <- cov_names
# hcdata <- dendro_data(hc)
#
# p_dendro <- ggplot() +
#   geom_segment(
#     data = segment(hcdata),
#     aes(x = x, y = y, xend = xend, yend = yend)
#   ) +
#   scale_x_continuous(
#     position = "top", breaks = hcdata$labels$x, labels = hcdata$labels$label
#   ) +
#   scale_y_reverse() +
#   ylab("$dist{theta}{theta^prime}$") +
#   xlab("Covariate ($theta$)") +
#   coord_flip() +
#   theme_bw() +
#   theme(panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.grid.minor.x = element_blank()
#   )
