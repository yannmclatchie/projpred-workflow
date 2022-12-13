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

get_case_study_weakly_relevant_sh_file <- function(time,R2,rho,fix_beta){
  experiment_name <- paste0('varsel_objects_R2_',R2,'_rho_',rho,'_fix_beta',fix_beta)
  sh_path <- paste0(experiment_name,'.sh')
  sh_file_str <- create_sh_file(mem='10G',time=time,output=paste0(experiment_name,'.out'),
                                script=paste0('R/case_study_weakly_relevant.R',
                                              ' --R2 ',R2,
                                              ' --rho ',rho,
                                              ' --fix_beta ',fix_beta))
  cat(sh_file_str,file=sh_path)
}

R2_vec <- c(0.3,0.4,0.5,0.7)
rho_vev <- c(0,0.1)
fix_beta_vec <- c(FALSE,TRUE)
batch_time <- '50:00:00'
count=0

for(R2 in R2_vec){
  for(rho in rho_vev){
    for(fix_beta in fix_beta_vec){
      count=count+1
      get_case_study_weakly_relevant_sh_file(time=batch_time,R2=R2,rho=rho,fix_beta=fix_beta)
    }
  }
}
