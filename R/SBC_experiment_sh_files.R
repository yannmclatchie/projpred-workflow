create_sh_file <- function(mem,time,output,script){
  str <- paste0('#!/bin/bash -l\n',
                '#SBATCH --time=',time,'\n',
                '#SBATCH --mem=',mem,'\n',
                '#SBATCH --cpus-per-task=4\n',
                '#SBATCH --output=',output,'\n',
                '\n',
                'module load r/4.1.1-python3\n',
                'srun Rscript ',script)
  return(str)
}

get_SBC_experiment_sh_file <- function(time,N_sim,N,n_rel,n_irrel,prior_ref,validate_search,experiment_suffix,path){
  experiment_name <- paste0(N_sim,'sim_',
                            N,'N_',
                            n_rel,'rel_',
                            n_irrel,'irrel_',
                            prior_ref,'prior',
                            ifelse(experiment_suffix!='',paste0('_',experiment_suffix),''))
  sh_dir <- 'R/SBC_sh_files/'
  if(!dir.exists(file.path(sh_dir))){
    dir.create(file.path(sh_dir))
  }
  sh_path <- paste0(sh_dir,experiment_name,'.sh')
  sh_file_str <- create_sh_file(mem='5G',time=time,output=paste0(experiment_name,'.out'),
                                script=paste0('../run_SBC_experiment.R',
                                              ' --N_sim ',N_sim,
                                              ' --N ',N,
                                              ' --n_rel ',n_rel,' --n_irrel ',n_irrel,
                                              ' --prior_ref ',prior_ref,
                                              ifelse(experiment_suffix!='',paste0(' --suffix ',experiment_suffix),''),
                                              ' --path ',path))
  cat(sh_file_str,file=sh_path)
}

N_sim <- 300
N <- 100
n_rel <- 10
n_irrel <- 60
prior_refs <- c('normal','normalr2d2')
batch_time <- '24:00:00'
batch_size <- 5
suffix='batch'
count=0
path <- 'results/'

for(pri in prior_refs){
  for(i in 1:(N_sim/batch_size)){
    count=count+1
    get_SBC_experiment_sh_file(time=batch_time,N=N,n_rel=n_rel,n_irrel=n_irrel,prior_ref=pri,N_sim=batch_size,experiment_suffix=paste(suffix,i,sep='_'),path=path)
  }
}

