.libPaths('../R/x86_64-pc-linux-gnu-library/4.1/') #tell where to load libraries
library(brms)
library(projpred)
source('R/generate.R')
set.seed(4321)
rho <- 0.9
n <- 500
p <- 100
dat <- .generate(rho=rho,n=n)
prior_ref <- prior(R2D2(mean_R2=0.3,prec_R2=5,cons_D2=1))
mod_ref <- brm(as.formula(paste0('y~',paste(paste0('X.',1:p),collapse='+'))),data=dat,prior=prior_ref,family=gaussian(),refresh=0,iter=2000,cores=4)

# rhats <- rhat(mod_ref)
# loo_ref <- loo(mod_ref)
# neff_ratio(mod_ref)
# elpd_ref <- loo_ref$estimates["elpd_loo", "Estimate"]

vs_forward_loo <- cv_varsel(mod_ref,cv_method='LOO',method = "forward",validate_search=F,nterms_max = 30)
vs_L1_loo <- cv_varsel(mod_ref,cv_method='LOO',method = "L1",validate_search=F,nterms_max = 30)
vs_forward_kfold <- cv_varsel(mod_ref,cv_method='kfold',method = "forward",nterms_max = 30)
size_maximum_elpd <- vs_forward_loo$summary$size[which.max(vs_forward_loo$summary$elpd.loo)]
print(paste('running cross validation of search path with max terms =',size_maximum_elpd))
vs_forward_loo_validated <- cv_varsel(mod_ref,cv_method='LOO',method = "forward",validate_search=T,nterms_max = size_maximum_elpd)
save(vs_forward_loo,vs_L1_loo,vs_forward_kfold,vs_forward_loo_validated,file=paste0('R/dat/case_study_4.1_',n,'.RData'))
