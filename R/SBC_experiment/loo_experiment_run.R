library(brms)
library(dplyr)
options(brms.backend = "rstan")
source('R/SBC_experiment/SBC_experiment_functions.R')
run_loo_experiment <- function(N_sim,N,n_rel,n_irrel){
  covars <- paste0('x',seq_len(n_rel+n_irrel))
  rel_covars <- paste0('x',seq_len(n_rel))
  rel_covars_str <- paste(rel_covars,collapse='+')
  ref_mod_vars_str <- paste(c('t',covars),collapse='+')
  ref_mod_formula <- paste0('y~',ref_mod_vars_str)
  loo_results <- tibble()
  for(i in 1:N_sim){
    print(paste('Iteration:',i))
    seed_iter <- sample(1:1e8,size = 1)
    print(paste0('seed: ',seed_iter))
    sim <- def_indep_treatment_experiment(N = N,n_rel = n_rel, n_irrel = n_irrel, seed = seed_iter)
    dat <- sim$dat
    beta <- sim$beta
    sigma <- sim$sigma
    #compile reference model in first iteration
    if(i==1){
        #reference model with independent normal prior
        prior_mod_normal <- prior(normal(0,1),class="b") + prior(exponential(1),class='sigma')
        mod_ref_normal <- brm(ref_mod_formula,data=dat,prior=prior_mod_normal,family=gaussian(),refresh=0,iter=2000,cores=4, save_pars = save_pars(all = TRUE))
        # rstan::check_hmc_diagnostics(mod_ref_normal$fit)
        #reference model with r2d2 prior on all coefficients
        #prior_mod_r2d2 <- prior(R2D2(mean_R2=0.3,prec_R2=20,cons_D2 = 0.2)) + prior(exponential(1),class='sigma')
        prior_mod_r2d2 <- prior(R2D2(mean_R2=0.3,prec_R2=20,cons_D2 = 0.2)) + prior(exponential(1),class='sigma')
        mod_ref_r2d2 <- brm(ref_mod_formula,data=dat,prior=prior_mod_r2d2,family=gaussian(),refresh=0,chains=4,iter=4000,thin=4,cores=4,control=list(adapt_delta=0.999), save_pars = save_pars(all = TRUE))
        # rstan::check_hmc_diagnostics(mod_ref_r2d2$fit)

        #reference model with normal prior on treatment and r2d2 prior on the rest
        mod_ref_r2d2normal<- brm(ref_mod_formula,data=dat,prior=prior_mod_r2d2,family=gaussian(),refresh=0,chains=0,iter=4000,thin=4,cores=4,control=list(adapt_delta=0.999), save_pars = save_pars(all = TRUE))
        mod_ref_standata <- make_standata(ref_mod_formula,data=dat,prior=prior_mod_r2d2)
        mod_ref_modified <- rstan::stan(file = 'stan/r2d2_normal.stan', data = mod_ref_standata,refresh=0)
        # rstan::check_hmc_diagnostics(mod_ref_modified)
        mod_ref_r2d2normal$fit <- mod_ref_modified
        mod_ref_r2d2normal <- rename_pars(mod_ref_r2d2normal) # note: warning message is because one concentration parameter was removed. Model object stays intact
        #store the compiled models to be used later
        mod_ref_normal_base <- mod_ref_normal
        mod_ref_r2d2_base <- mod_ref_r2d2
        mod_ref_r2d2normal_base <- mod_ref_r2d2normal
    }else{
      mod_ref_normal <- update(mod_ref_normal_base,newdata=dat,recompile=F,refresh=0)
      mod_ref_r2d2 <- update(mod_ref_r2d2_base,newdata=dat,recompile=F,refresh=0)
      mod_ref_r2d2normal <- update(mod_ref_r2d2normal_base,newdata=dat,recompile=F,refresh=0)
    }
    loo_mod_ref_normal <- loo(mod_ref_normal, moment_match = TRUE, reloo = TRUE)
    loo_mod_ref_r2d2 <- loo(mod_ref_r2d2, reloo = TRUE)
    loo_mod_ref_r2d2normal <- loo(mod_ref_r2d2normal, reloo = TRUE)
    loo_normal_vs_r2d2 <- loo_compare(loo_mod_ref_normal,loo_mod_ref_r2d2)
    loo_normal_vs_r2d2normal <- loo_compare(loo_mod_ref_normal,loo_mod_ref_r2d2normal)
    r2=summary(lm(as.formula(paste0('y~t+',paste(paste0('x',1:n_rel),collapse='+'))),data=dat))$r.squared
    loo_results <- bind_rows(loo_results,tibble(iter=i,
                                         r2=r2,
                                         elpd_ref_normal=loo_mod_ref_normal$estimates['elpd_loo','Estimate'],
                                         normal_r2d2_loser=rownames(loo_normal_vs_r2d2)[2],
                                         elpd_diff_normal_r2d2=loo_normal_vs_r2d2[2,'elpd_diff'],
                                         se_diff_normal_r2d2=loo_normal_vs_r2d2[2,'se_diff'],
                                         normal_r2d2normal_loser=rownames(loo_normal_vs_r2d2normal)[2],
                                         elpd_diff_normal_r2d2normal=loo_normal_vs_r2d2normal[2,'elpd_diff'],
                                         se_diff_normal_r2d2normal=loo_normal_vs_r2d2normal[2,'se_diff']))
  }
  return(loo_results)
}


#Try running one iteration of loo experiment
set.seed(1234)
loo_results <- run_loo_experiment(N_sim=10,N=100,n_rel=10,n_irrel=60)
write.table(loo_results,'loo_experiment_r2d2_0.3_20_0.2.tsv',quote = F,row.names=F,sep = '\t')
