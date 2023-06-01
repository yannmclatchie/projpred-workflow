## ----- Helper functions -------
def_indep_treatment_experiment <- function(N,n_rel,n_irrel,seed,sd_beta=1,lambda=1){
  set.seed(seed)
  n_covars <- n_rel + n_irrel
  covars <- paste0('x',seq_len(n_covars))
  beta <- rnorm(n_rel+1,mean=0,sd=sd_beta)
  names(beta) <- c('t',covars[seq_len(n_rel)])
  sigma <- rexp(n=1,rate=lambda)
  names(sigma) <- 'sigma'
  dat <- as.data.frame(matrix(rnorm(n_covars*N,mean=0,sd=1),nrow=N))
  names(dat) <- covars
  dat <- cbind(data.frame(t=rbinom(N,size=1,prob=0.5)),dat) # randomized trial
  dat$mean <- as.matrix(dat[,seq_len(n_rel+1)]) %*%beta
  dat$y <- rnorm(N,mean=dat$mean,sd=sigma)
  return(list(dat=dat,beta=beta,sigma=sigma))
}

get_posterior_summary <- function(post_dat){
  group_by(post_dat,sim,model) %>%
    summarise(prior_draw=unique(prior_draw),
              pit=mean(value < prior_draw),
              post_mean=mean(value),
              post_sd=sd(value),
              lower_95=quantile(value,0.025),
              upper_95=quantile(value,0.975),
              interval_width=upper_95-lower_95,
              squared_error=(prior_draw-post_mean)^2,
              cover=prior_draw<=upper_95 & prior_draw>=lower_95)
}


run_SBC_experiment <- function(N_sim,N,n_rel,n_irrel,prior_ref){
  covars <- paste0('x',seq_len(n_rel+n_irrel))
  rel_covars <- paste0('x',seq_len(n_rel))
  rel_covars_str <- paste(rel_covars,collapse='+')
  ref_mod_vars_str <- paste(c('t',covars),collapse='+')
  ref_mod_formula <- paste0('y~',ref_mod_vars_str)
  treatment_post_list <- list()
  sigma_post_list <- list()
  run_info_list <- list()
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
      if(prior_ref=='normal'){
        prior_mod <- prior(normal(0,1),class="b") + prior(exponential(1),class='sigma')
        mod_ref <- brm(ref_mod_formula,data=dat,prior=prior_mod,family=gaussian(),refresh=0,iter=2000,cores=4) # save_pars = save_pars(all = TRUE)
      }else if(prior_ref=='normalr2d2'){
        prior_mod <- prior(R2D2(mean_R2=0.3,prec_R2=20,cons_D2 = 0.2)) + prior(exponential(1),class='sigma')
        mod_ref <- brm(ref_mod_formula,data=dat,prior=prior_mod,family=gaussian(),refresh=0,chains=0,iter=16000,thin=8,cores=4,control=list(adapt_delta=0.999)) # save_pars = save_pars(all = TRUE)
        mod_ref_standata <- make_standata(ref_mod_formula,data=dat,prior=prior_mod)
        mod_ref_modified <- rstan::stan(file = 'stan/r2d2_normal.stan', data = mod_ref_standata,refresh=0)
        mod_ref$fit <- mod_ref_modified
        mod_ref <- rename_pars(mod_ref)
      }else{
        stop('Prior not supported')
      }
      mod_ref_base <- mod_ref
      #fit model with intercept and treatment only. Is needed in case variable selection chooses treatment only and prior is R2D2
      prior_trivial <- prior(normal(0,1),class='b') + prior(exponential(1),class='sigma')
      mod_trivial <- brm(y~t,data=dat,prior=prior_trivial,family=gaussian(),refresh=0,iter=2000,cores=4)
      mod_trivial_base <- mod_trivial
    }else{
      mod_ref <- update(mod_ref_base,newdata=dat,recompile=F,refresh=0)
    }
    #calculate loo diagnostics for reference model
    loo_mod_ref <- loo(mod_ref)
    ### For actual lpds:
    # loo_mod_ref_clean <- loo(mod_ref, moment_match = TRUE, reloo = TRUE)
    ###
    # TODO: Does the current code allow to compute elpd difference between different `prior_ref`s but for the same dataset, i.e., the same SBC iteration (after storing `loo_mod_ref_clean` in the results, I guess)?
    #Projpred variable selection using kfold cross validation
    cv_select_prj <- cv_varsel(mod_ref,cv_method='kfold',method='forward',validate_search=T,nterms_max=30,search_terms=paste('t + ',covars))
    selected_vars_prj <- sort(solution_terms(cv_select_prj)[seq_len(suggest_size(cv_select_prj))])
    #In case nothing was selected, add t to enforce choosing t in the model
    selected_vars_prj <- unique(c('t',unlist(strsplit(selected_vars_prj,split=' \\+ '))))
    selected_vars_prj_str <- paste(selected_vars_prj,collapse='+')
    #project on selected variables from cross validation
    prj_sel <- project(mod_ref, solution_terms = selected_vars_prj, ndraws=4000)
    #Model refit
    if(selected_vars_prj_str=='t'){
      mod_refit_prj <- update(mod_trivial_base,newdata=dat,recompile=F,refresh=0)
    }else{
      mod_refit_prj <- update(mod_ref_base,formula.=as.formula(paste0('y~',selected_vars_prj_str,collapse='+')),newdata=dat,recompile=F,refresh=0)
    }
    mod_names <- c('mod_ref','prj_sel','mod_refit_prj')
    mods <- list(mod_ref,prj_sel,mod_refit_prj)

    treatment_post_list[[i]] <- lapply(mods,function(m) tibble(value=as.matrix(m)[,'b_t'])) %>%
                         bind_rows(.id='model_nr') %>%
                         mutate(model=mod_names[as.integer(model_nr)],
                                sim=i,
                                prior_draw=beta['t']) %>%
                         dplyr::select(sim,model,prior_draw,value)
    sigma_post_list[[i]] <- lapply(mods,function(m) tibble(value=as.matrix(m)[,'sigma'])) %>%
                            bind_rows(.id='model_nr') %>%
                            mutate(model=mod_names[as.integer(model_nr)],
                                   sim=i,
                                   prior_draw=sigma) %>%
                            dplyr::select(sim,model,prior_draw,value)

    #Store model diagnostics of the reference model as well as other run info
    mod_ref_summary <- summary(mod_ref)
    run_info_list[[i]] <- bind_cols(tibble(sim=i,seed=seed_iter),
                                    tibble(snr=sigma^2/var(dat$y),
                                           max_r_hat = min(mod_ref_summary$fixed$Rhat),
                                           min_bulk_ess_ref = min(mod_ref_summary$fixed$Bulk_ESS),
                                           min_tail_ess_ref = min(mod_ref_summary$fixed$Tail_ESS),
                                           num_divergent_ref = rstan::get_num_divergent(mod_ref$fit),
                                           num_high_k_value = sum(loo_mod_ref$diagnostics$pareto_k>0.7)),
                                    as_tibble(t(c(beta,sigma))),
                                    tibble(full_mod=ref_mod_vars_str,
                                           selected_mod_prj=selected_vars_prj_str))
  }
  treatment_post_dat <- bind_rows(treatment_post_list)
  sigma_post_dat <- bind_rows(sigma_post_list)
  treatment_summary <- get_posterior_summary(treatment_post_dat)
  sigma_summary <- get_posterior_summary(sigma_post_dat)
  run_info <- bind_rows(run_info_list)
  return(list(treatment_summary=treatment_summary,sigma_summary=sigma_summary,run_info=run_info))
}
