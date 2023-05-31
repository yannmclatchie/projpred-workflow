library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(vroom)
library(bayesplot)
library(gridExtra)
library(grid)
library(patchwork)
source('R/aux/aux_plotting.R')

plot_calibration <- function(dat){
  ppc_pit_ecdf_grouped(
    pit = dat$pit,
    group = dat$var,
    interpolate_adj = F) +
    theme(
      plot.title=element_text(size=12,family=''),
      strip.text.x = element_blank()
    )
}

create_layout_row <- function(row_id,plot_ids,width,height){
  plot_ids_rep <- paste(sapply(plot_ids,function(id) rep(id,width)),collapse='')
  return(paste(sapply(1:height,function(i) paste0(row_id,plot_ids_rep)),collapse='\n'))
}

pri <- c('normal','normalr2d2')
mod_names <- c('Reference model'='mod_ref','Projected model'='prj_sel')
cols <- c('#BC3C29FF','#0072B5FF','#E18727FF','#20854EFF','#7876B1FF','#6F99ADFF')

#read in the result files from the SBC experiment
results_path <- 'data/SBC_experiment/'
results_files <- list.files(results_path,full.names=T)
treatment_summary <- lapply(results_files[grepl('treatment_summary',results_files)],function(f){
  prior=ifelse(grepl('r2d2',f),'normalr2d2','normal')
  read_tsv(f) %>% mutate(prior=prior,batch=as.numeric(gsub('\\.tsv','',gsub('^.*batch_','',f))))
}) %>% bind_rows()

run_info <- lapply(results_files[grepl('run_info',results_files)],function(f){
  prior=ifelse(grepl('r2d2',f),'normalr2d2','normal')
  read_tsv(paste0(f)) %>% mutate(prior=prior,batch=as.numeric(gsub('\\.tsv','',gsub('^.*batch_','',f))))
}) %>% bind_rows()
#select iterations where the posterior sampling was of sufficient quality
run_info_filtered <- filter(run_info,min_bulk_ess_ref>400,min_tail_ess_ref>400,num_divergent_ref<=5)
set.seed(12)
treatment_summary_filtered <- inner_join(treatment_summary,
                                         select(run_info_filtered,prior,sim,batch,min_bulk_ess_ref,min_tail_ess_ref,num_divergent_ref),
                                         by=c('prior','batch','sim')) %>%
                              group_by(batch,sim) %>%
                              mutate(group_nr=cur_group_id()) %>%
                              ungroup() %>%
                              group_by(prior) %>%
                              filter(group_nr %in% sample(unique(group_nr),size=300,replace=F))


plot_list <- lapply(pri,function(p){
  pri_dat <- filter(treatment_summary_filtered,prior==p)
  lapply(1:length(mod_names),function(i){
    pri_mod_dat <- filter(pri_dat,model==mod_names[i])
    ppc_pit_ecdf(pit = pri_mod_dat$pit,interpolate_adj = F,plot_diff=T) +
    ggtitle(ifelse(p=='normal',names(mod_names)[i],'')) +
    theme_classic()
  })
}) %>% unlist(recursive = F)

rows_titles <- lapply(c('Gaussian prior','R2D2 prior'),function(m){
  ggplot() + annotate(geom = 'text', x=1, y=1, label=m, angle = 90,size=4.5,family= theme_get()$text[["family"]]) + theme_void()
})

ids <- letters[seq_len(length(rows_titles)+length(plot_list))]
layout_lines <- sapply(1:length(rows_titles),function(i){
  id_idx <- seq((length(rows_titles)+1+(i-1)*length(mod_names)),length(rows_titles)+i*length(mod_names ))
  create_layout_row(row_id=ids[i],plot_ids=ids[id_idx],height=4,width=4)
})
layout <- paste0('\n',paste(layout_lines,collapse='\n'),'\n')

patches <- c(rows_titles,plot_list)
names(patches) <- ids

sbc_panel <- wrap_plots(patches, guides = 'collect', design = layout)
save_tikz_plot(sbc_panel, width = 6, filename = "tex/calibration_plot.tex")
