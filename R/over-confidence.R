library(ggplot2)
library(dplyr)
library(tidyr)
library(brms)
library(projpred)
library(patchwork)
source('R/aux/aux_plotting.R')

#### ---- Define experiment and simulate data ---- 
seed <- 6089090
set.seed(seed)
N <- 100
n_rel <- 10
n_irrel <- 80
n_covars <- n_rel + n_irrel
covars <- paste0('x',1:n_covars)
beta <- rnorm(n_rel,mean=0,sd=1) 
sigma <- rexp(1)
dat <- as.data.frame(matrix(rnorm(n_covars*N,mean=0,sd=1),nrow=N))
names(dat) <- covars
dat$y <- rnorm(N,mean=as.matrix(dat[,1:n_rel]) %*%beta,sd=sigma)
ref_mod_vars_str <- paste(covars,collapse='+') 
prior_mod <- prior(normal(0,1),class="b") + prior(exponential(1),class='sigma')
mod_ref <- brm(paste0('y~',ref_mod_vars_str),data=dat,prior=prior_mod,family=gaussian(),refresh=0,iter=4000,cores=4)

DGP_covars <- covars[1:n_rel]
DGP_covars_str <- paste(DGP_covars,collapse='+')
large_covars <- c(DGP_covars,paste0('x',seq(n_rel+1,n_rel+10)))
large_covars_str <- paste(large_covars,collapse='+')
small_covars <- c(paste0('x',seq(1,floor(n_rel/2))))
small_covars_str <- paste(small_covars,collapse='+')

prj_large <- project(mod_ref,solution_terms=large_covars,ndraws=8000)
prj_DGP <- project(mod_ref,solution_terms=DGP_covars,ndraws=8000)
prj_small <- project(mod_ref,solution_terms=small_covars,ndraws=8000)

mod_refit_large <- brm(as.formula(paste0('y ~ ',large_covars_str)),data=dat,prior=prior_mod,family=gaussian(),refresh=0,iter=4000,cores=4)
mod_refit_DGP <- brm(as.formula(paste0('y ~ ',DGP_covars_str)),data=dat,prior=prior_mod,family=gaussian(),refresh=0,iter=4000,cores=4)
mod_refit_small <- brm(as.formula(paste0('y ~ ',small_covars_str)),data=dat,prior=prior_mod,family=gaussian(),refresh=0,iter=4000,cores=4)

#### ---- Process results ----
extract_post <- function(mod,mod_name){
    post_mat <- as.matrix(mod)[,c('b_x1','b_x2','sigma')]
    colnames(post_mat) <- paste0(mod_name,':',colnames(post_mat))
    as_tibble(post_mat)
}
mod_types <- c('prj'='Projection','refit'='Refit','ref'='Reference model')
var_names <- c('b_x1'='beta[1]','b_x2'='beta[2]','sigma'='sigma')
mod_list <- list(mod_ref,prj_large,prj_DGP,prj_small,mod_refit_large,mod_refit_DGP,mod_refit_small)
posterior_dat <- bind_cols(extract_post(mod_ref,'ref_all'),
                           extract_post(prj_large,'prj_large'),
                           extract_post(prj_DGP,'prj_DGP'),
                           extract_post(prj_small,'prj_small'),
                           extract_post(mod_refit_large,'refit_large'),
                           extract_post(mod_refit_DGP,'refit_DGP'),
                           extract_post(mod_refit_small,'refit_small'),) %>%
    pivot_longer(cols=everything(),names_to='var',values_to='value') %>%
    separate(var,sep = ':',into = c('model','var'),remove = F) %>%
    separate(model,sep = '_',into = c('model_type','size'),remove = F) %>%
    mutate(mod_type_name=mod_types[model_type],
           var_name=var_names[var])


#### ---- Plot ----
cols <- c('#BC3C29FF','#0072B5FF','#E18727FF','#20854EFF','#7876B1FF','#6F99ADFF','#EE4C97FF')[c(3,1,2)]

true_vals <- tibble(var_name=var_names,value=c(beta[1:2],sigma))
large_mod_p <- ggplot(filter(posterior_dat,size %in% c('all','large'))) +
                geom_density(aes(value,col=model,fill=model),alpha=0.3) +
                facet_wrap(~var_name,scale='free',labeller=label_parsed) +
                geom_vline(data=true_vals,aes(xintercept=value),linetype='dashed') +
                scale_y_continuous(expand=c(0,0)) +
                scale_color_manual(values=cols,name='Larger model') +
                scale_fill_manual(values=cols,name='Larger model') +
                xlab('') +
                ylab('Density') +
                theme_classic() +
                theme(legend.position=c(0.9,0.8),
                      strip.background=element_blank(),
                      strip.text=element_text(size=12))
ggsave('fig/cmp_larger_model.png',plot=large_mod_p,width=8,height=4)

DGP_mod_p <- ggplot(filter(posterior_dat,size %in% c('all','DGP'))) +
                geom_density(aes(value,col=model,fill=model),alpha=0.3) +
                facet_wrap(~var_name,scale='free',labeller=label_parsed) +
                geom_vline(data=true_vals,aes(xintercept=value),linetype='dashed') +
                scale_y_continuous(expand=c(0,0)) +
                scale_color_manual(values=cols,name='DGP') +
                scale_fill_manual(values=cols,name='DGP') +
                xlab('') +
                ylab('Density') +
                theme_classic() +
                theme(legend.position=c(0.9,0.8),
                      strip.background=element_blank(),
                      strip.text=element_text(size=12))
ggsave('fig/cmp_DGP_model.png',plot=DGP_mod_p,width=8,height=4)

small_mod_p <- ggplot(filter(posterior_dat,size %in% c('all','small'))) +
                geom_density(aes(value,col=model,fill=model),alpha=0.3) +
                facet_wrap(~var_name,scale='free',labeller=label_parsed) +
                geom_vline(data=true_vals,aes(xintercept=value),linetype='dashed') +
                scale_y_continuous(expand=c(0,0)) +
                scale_color_manual(values=cols,name='Smaller model') +
                scale_fill_manual(values=cols,name='Smaller model') +
                xlab('') +
                ylab('Density') +
                theme_classic() +
                theme(legend.position=c(0.9,0.8),
                      strip.background=element_blank(),
                      strip.text=element_text(size=12))
ggsave('fig/cmp_small_model.png',plot=small_mod_p,width=8,height=4)

p1 <- DGP_mod_p +
      ggtitle('Oracle model') +
      theme(legend.position='none',
            plot.title=element_text(size=12))

p2 <- small_mod_p +
      ggtitle('Smaller model') +
      theme(legend.position='none',
            plot.title=element_text(size=12))

combined_p <- p1 / p2

ggsave('fig/overconfidence.png',plot=combined_p,width=8,height=5)

save_tikz_plot(combined_p,width=6,height=4,filename = 'tex/overconfidence.tex')


