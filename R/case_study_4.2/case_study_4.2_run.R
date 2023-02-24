.libPaths('../R/x86_64-pc-linux-gnu-library/4.1/') #tell which libraries to use
library(brms)
library(projpred)
library(argparser)
source('R/generate.R')
## ----- Define argument parser -------
pars <- arg_parser("weakly relevant predictors case study from the paper 'Robust and efficient projection predictive inference'")
pars <- add_argument(pars, "--R2", help="R2 in the simulated data",type="numeric")
pars <- add_argument(pars, "--fix_beta", help="Chould betas all be equal?",type="logical")
pars <- add_argument(pars, "--rho", help="correlation coefficients between the predictors",type="numeric")

args <- parse_args(pars)
set.seed(4321)
t <- Sys.time()
rho <- args$rho
n <- 500
p <- 50
R2 <- args$R2
sigma <- 1
sim <- generate_fixed_R2_data(n=n,p=p,sigma=sigma,R2=args$R2,rho=args$rho,fixed_beta = args$fix_beta)
dat <- sim$df
beta <- sim$beta
prior_ref <- prior(R2D2(mean_R2=0.3,prec_R2=5,cons_D2=10))
mod_ref <- brm(as.formula(paste0('y~',paste(paste0('X',1:p),collapse='+'))),data=dat,prior=prior_ref,family=gaussian(),refresh=0,iter=2000,cores=4)

# rhats <- rhat(mod_ref)
# loo_ref <- loo(mod_ref)
# neff_ratio(mod_ref)
# elpd_ref <- loo_ref$estimates["elpd_loo", "Estimate"]

vs_forward_loo <- cv_varsel(mod_ref,cv_method='LOO',method = "forward",validate_search=F,nterms_max = p)
size_maximum_elpd <- vs_forward_loo$summary$size[which.max(vs_forward_loo$summary$elpd.loo)]
print(paste('Running cross validated search with nterms_max = ',size_maximum_elpd))
vs_forward_loo_validated <- cv_varsel(mod_ref,cv_method='LOO',method = "forward",validate_search=T,nterms_max = size_maximum_elpd)
save(beta,vs_forward_loo,vs_forward_loo_validated,
     file=paste0('results/case_study_4.2_R2_',args$R2,'_rho_',args$rho,'_fix_beta',args$fix_beta,'.RData'))
print(Sys.time()-t)
