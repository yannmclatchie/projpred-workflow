.libPaths('../R/x86_64-pc-linux-gnu-library/4.1/')
library(brms)
library(projpred)
library(argparser)
source('R/generate.R')
#source('R/aux/aux_plotting.R')
## ----- Define argument parser -------
p <- arg_parser("weakly relevant predictors case study from the paper 'Robust and efficient projection predictive inference'")
p <- add_argument(p, "--R2", help="R2 in the simulated data",type="numeric")
p <- add_argument(p, "--fix_beta", help="Chould betas all be equal?",type="logical")
p <- add_argument(p, "--rho", help="correlation coefficients between the predictors",type="numeric")

args <- parse_args(p)
set.seed(4321)
rho <- args$rho
n <- 500
p <- 100
R2 <- args$R2
sigma <- 1
sim <- generate_fixed_R2_data(n=n,p=p,sigma=sigma,R2=args$R2,rho=args$rho,fixed_beta = args$fix_beta)
dat <- sim$df
beta <- sim$beta
prior_ref <- prior(R2D2(mean_R2=0.3,prec_R2=5,cons_D2=10))
mod_ref <- brm(as.formula(paste0('y~',paste(paste0('X',1:100),collapse='+'))),data=dat,prior=prior_ref,family=gaussian(),refresh=0,iter=2000,cores=4)

# rhats <- rhat(mod_ref)
# loo_ref <- loo(mod_ref)
# neff_ratio(mod_ref)
# elpd_ref <- loo_ref$estimates["elpd_loo", "Estimate"]

vs_forward_loo <- cv_varsel(mod_ref,cv_method='LOO',method = "forward",validate_search=F,nterms_max = 100)
size_maximum_elpd <- vs_forward_loo$summary$size[which.max(vs_forward_loo$summary$elpd.loo)]
vs_forward_loo_validated <- cv_varsel(mod_ref,cv_method='LOO',method = "forward",validate_search=T,nterms_max = size_maximum_elpd)
save(vs_forward_loo,vs_forward_loo_validated,
     file=paste0('results/varsel_objects_R2_',args$R2,'_rho_',args$rho,'_fix_beta',args$fix_beta,'.RData'))

# mod_size_candidates <- c('confidence_heuristic'=suggest_size(vs_forward_loo),
#                         'elpd_difference_heuristic'=suggest_size(vs_forward_loo,thres_elpd = -4))
#
# vs_forward_loo$summary <- mutate(vs_forward_loo$summary,mod_size_candidates=size %in% mod_size_candidates)
#
# ### ELPD plot
# ggplot(data=vs_forward_loo$summary,
#        aes(x=size,
#            y=elpd.loo,
#            ymin=elpd.loo-se,
#            ymax=elpd.loo+se,
#            col=mod_size_candidates)) +
# geom_pointrange() +
# geom_line() +
# geom_hline(yintercept=elpd_ref,colour = "red", linetype = "longdash") +
# geom_vline(xintercept = mod_size_candidates['confidence_heuristic'], colour = "grey", linetype = "longdash") +
# geom_vline(xintercept = mod_size_candidates['elpd_difference_heuristic'], colour = "grey", linetype = "longdash") +
# annotate("text", x = 1.5, y = elpd_ref + 10, colour = "red", label = "Reference model elpd",size=3) +
# scale_color_manual(values=c('blue','black'),name='') +
# scale_x_continuous(breaks = 0:max_size_displayed) +
# xlab('Model size') +
# ylab('$elpd$') +
# theme_bw() +
# theme(panel.grid.major.x = element_blank(),
#       panel.grid.minor.x = element_blank(),
#       panel.grid.minor.y = element_blank(),
#       legend.position = c(0.8,0.55))
# ggsave('R/fig/elpd_plot_case_study_3_1.png',elpd_plot,width=6,height=4)
# save_tikz_plot(elpd_plot,width=6,height=4,filename = 'tex/elpd_plot_case_study_3_1.tex')
