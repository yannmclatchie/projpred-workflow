## ----- Setup -------
#move from sbatch directory to project directory
setwd('../../')
#Tell cluster where R is located
.libPaths('../R/x86_64-pc-linux-gnu-library/4.1/')
library(brms)
library(projpred)
library(dplyr)
library(argparser)
source('R/SBC_experiment/SBC_experiment_functions.R')
source('R/SBC_experiment/loo_experiment_functions.R')
## ----- Define argument parser -------
p <- arg_parser("Run loo for the SBC experiment from the paper 'Robust and efficient projection predictive inference'")
p <- add_argument(p, "--N_sim", help="Number of SBC simulations to perform", default=10,type="numeric")
p <- add_argument(p, "--N", help="Integer determining the size of the simulated data sets",type="numeric")
p <- add_argument(p, "--n_rel", help="Integer determining the number of relevant predictors",type="numeric")
p <- add_argument(p, "--n_irrel", help="Integer determining the number of irrelevant predictors",type="numeric")
p <- add_argument(p, "--prior_ref", help="Character denoting prior of the beta coefficients in the reference model.
                                          Supports values 'normal' and 'normalr2d2', where normal r2d2 refers to a
                                          normal prior on treatment and r2d2 on remainin covariates", default='normal',type="character")
p <- add_argument(p, "--experiment_suffix", help="Suffix to add to file path to label experiment.", default="",type="character")
p <- add_argument(p, "--path", help="Path were results should be saved", default="",type="character")

args <- parse_args(p)

## ----- Run experiment -------
# args <- list()
# args$N_sim <- 1
# args$N <- 100
# args$n_rel <- 10
# args$n_irrel <- 10
# args$prior_ref <- 'normal'
# args$experiment_suffix <- ''
# args$path <- 'results/'

experiment <- run_loo_experiment(N_sim=args$N_sim,N=args$N,n_rel=args$n_rel,n_irrel=args$n_irrel,prior_ref=args$prior_ref)

experiment_name <- paste0(args$N_sim,'sim_',
                          args$N,'N_',
                          args$n_rel,'rel_',
                          args$n_irrel,'irrel_',
                          args$prior_ref,'prior',
                          ifelse(args$experiment_suffix!='',paste0('_',args$experiment_suffix),''))
write.table(experiment$treatment_summary,file=paste0(args$path,'treatment_summary_',experiment_name,'.tsv'),quote=F,row.names=F,sep='\t')
write.table(experiment$sigma_summary,file=paste0(args$path,'sigma_summary_',experiment_name,'.tsv'),quote=F,row.names=F,sep='\t')
write.table(experiment$run_info,file=paste0(args$path,'run_info_',experiment_name,'.tsv'),quote=F,row.names=F,sep='\t')
