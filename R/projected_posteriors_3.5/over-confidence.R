library(ggplot2)
library(dplyr)
library(tidyr)
library(brms)
library(projpred)
library(patchwork)
source('R/aux/aux_plotting.R')

# generate data with independent covariates, where target is a linear combination
# of n_rel relevant covariates, with an additional gaussian noise.
def_indep_experiment <- function(N,n_rel,n_irrel,seed){
    set.seed(seed)
    n_covars <- n_rel + n_irrel
    covars <- paste0('x',1:n_covars)
    beta <- rnorm(n_rel,mean=0,sd=1)
    sigma <- rexp(1)
    dat <- as.data.frame(matrix(rnorm(n_covars*N,mean=0,sd=1),nrow=N))
    names(dat) <- covars
    dat$y <- rnorm(N,mean=as.matrix(dat[,1:n_rel]) %*%beta,sd=sigma)
    return(list(dat=dat,covars=covars,beta=beta,sigma=sigma))
}

#### ---- Define experiment and simulate data ----
seed <- 6089090
N <- 100
n_rel <- 15
n_irrel <- 80
sim <- def_indep_experiment(N=N,n_rel=n_rel,n_irrel=n_irrel,seed=seed)
dat <- sim$dat
covars <- sim$covars
beta <- sim$beta
sigma <- sim$sigma
ref_mod_vars_str <- paste(covars,collapse='+')
prior_mod <- prior(normal(0,1),class="b") + prior(exponential(1),class='sigma')
mod_ref <- brm(paste0('y~',ref_mod_vars_str),data=dat,prior=prior_mod,family=gaussian(),refresh=0,iter=4000,cores=4)
SNR <- var(as.matrix(dat[,1:n_rel]) %*%beta)/sigma^2

irrel_cor <- sapply(seq(n_rel+1,n_rel+n_irrel),function(i) cor(dat$y,dat[,i]))
#covariates of a model with fewer covariates, with half of the relevant variables
small_covars <- c(paste0('x',seq(1,n_rel)))#c(paste0('x',seq(1,floor(n_rel/2))))
small_covars_str <- paste(small_covars,collapse='+')

prj_small <- project(mod_ref,solution_terms=small_covars,ndraws=8000)
mod_refit_small <- brm(as.formula(paste0('y ~ ',small_covars_str)),data=dat,prior=prior_mod,family=gaussian(),refresh=0,iter=4000,cores=4)

#### ---- Process results ----
extract_post <- function(mod,mod_name){
    post_mat <- as.matrix(mod)[,c('b_x1','b_x2','sigma')]
    colnames(post_mat) <- paste0(mod_name,':',colnames(post_mat))
    as_tibble(post_mat)
}

extract_post_pred <- function(mod,mod_type){
    if(inherits(mod,'brmsfit')){
        pp <- posterior_predict(mod)
    }else if(inherits(mod,'projection')){
        pp <- proj_predict(mod)
    }else{
        stop("model type not recognized.")
    }
    as_tibble(pp) %>%
    mutate(mod_type_name=mod_type,
           sim=1:n()) %>%
    pivot_longer(cols=-c(sim,mod_type_name),names_to='datapoint',values_to='value')
}
mod_types <- c('ref'='Reference model','refit'='Submodel refit','prj'='Projected model')
var_names <- c('b_x1'='$beta_1$','b_x2'='$beta_2$','sigma'='$sigma$')
mod_list <- list(mod_ref,prj_small,mod_refit_small)
posterior_dat <- bind_cols(extract_post(prj_small,'prj_small'),
                           extract_post(mod_ref,'ref_all'),
                           #extract_post(mod_refit_small,'refit_small')
                           ) %>%
                pivot_longer(cols=everything(),names_to='var',values_to='value') %>%
                separate(var,sep = ':',into = c('model','var'),remove = F) %>%
                separate(model,sep = '_',into = c('model_type','size'),remove = F) %>%
                mutate(mod_type_name=mod_types[model_type],
                       var_name=var_names[var])

pp_summary <- bind_rows(extract_post_pred(prj_small,mod_types['prj']),
                        #extract_post_pred(mod_refit_small,mod_types['refit']),
                        extract_post_pred(mod_ref,mod_types['ref'])) %>%
                group_by(mod_type_name,datapoint) %>%
                summarise(lower=quantile(value,0.025),
                          mean = mean(value),
                          upper=quantile(value,0.975)) %>%
                mutate(var_name='$y$')

#### ---- Plot ----
cols <- c('#BC3C29FF','#0072B5FF','#E18727FF','#20854EFF','#7876B1FF','#6F99ADFF','#EE4C97FF')[seq_len(3)]
cols <- c("black", "grey", "blue")
true_vals <- tibble(var_name=var_names,value=c(beta[1:2],sigma))
small_post_plot_dat <- bind_rows(select(posterior_dat,mod_type_name,var_name,value),
                               select(pp_summary,mod_type_name,var_name,value=mean))
small_post_plot <- ggplot(small_post_plot_dat) +
  geom_density(data=mutate(dat,var_name='$y$'),aes(y),color=NA,fill='black',alpha=0.2,adjust=2) +
  stat_density(aes(value,col=mod_type_name),geom='line',position = "identity",adjust=2,size=0.8) +
  facet_wrap(~var_name,nrow=1,scale='free') +
  geom_vline(data=true_vals,aes(xintercept=value),linetype='dashed',size=0.7,alpha=1) +
  scale_y_continuous(expand=c(0.01,0)) +
  scale_color_manual(values=cols,breaks = mod_types,name='') +
  scale_fill_manual(values='black',breaks = 'Data density',name='') +
  xlab('') +
  ylab('Density') +
  theme_classic() +
  theme(legend.position='bottom',
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing = unit(0, "pt"),
        legend.text = element_text(size=12),
        strip.background=element_blank(),
        strip.text=element_text(size=12))
small_post_plot

ggsave('fig/overconfidence.png',plot=small_post_plot,width=8,height=3)
save_tikz_plot(small_post_plot,width=7,height=2.5,filename = 'tex/overconfidence.tex')
