library(ggplot2)
library(dplyr)
library(tidyr)
library(vroom)
library(bayesplot) # Note using version from github with SHA 9f8ff05df44faa57e4187fbebbefcdea68ae83cc
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
      strip.text.x = element_blank(),
      #plot.subtitle = element_text(size=16),
      #axis.title.x=element_text(size=18),
      #axis.title.y=element_text(size=18),
      #axis.text.x=element_text(size=16),
      #axis.text.y=element_text(size=16)
    )
}

create_layout_row <- function(row_id,plot_ids,width,height){
  plot_ids_rep <- paste(sapply(plot_ids,function(id) rep(id,width)),collapse='')
  return(paste(sapply(1:height,function(i) paste0(row_id,plot_ids_rep)),collapse='\n'))
}

pri <- c('normal','normalr2d2')
mod_names <- c('Reference model'='mod_ref','Submodel refit'='mod_refit_prj','Projedcted model'='prj_sel')
cols <- c('#BC3C29FF','#0072B5FF','#E18727FF','#20854EFF','#7876B1FF','#6F99ADFF')

treatment_summary <- read.table('R/dat/SBC_experiment_treatment_summary.tsv',header=T)
run_info <- read.table('R/dat/SBC_experiment_run_info.tsv',header=T)
run_info_filtered <- filter(run_info,min_bulk_ess_ref>400,min_tail_ess_ref>400,num_divergent_ref<=10)
treatment_summary_filtered <- inner_join(treatment_summary,
                                         select(run_info_filtered,prior,sim,batch,min_bulk_ess_ref,min_tail_ess_ref,num_divergent_ref),
                                         by=c('prior','batch','sim'))

plot_list <- lapply(pri,function(p){
  pri_dat <- filter(treatment_summary_filtered,prior==p)
  lapply(1:length(mod_names),function(i){
    pri_mod_dat <- filter(pri_dat,model==mod_names[i])
    plot_calibration(pri_mod_dat) +
      ggtitle(ifelse(p=='normal',names(mod_names)[i],''))
  })
}) %>% unlist(recursive = F)

rows_titles <- lapply(c('Over-fitting','Well-specified'),function(m){
  ggplot() + annotate(geom = 'text', x=1, y=1, label=m, angle = 90,size=4.5,family= theme_get()$text[["family"]]) + theme_void()
})

ids <- letters[seq_len(length(rows_titles)+length(plot_list))]
layout_lines <- sapply(1:length(rows_titles),function(i){
  id_idx <- seq((length(rows_titles)+1+(i-1)*3),length(rows_titles)+i*3)
  create_layout_row(row_id=ids[i],plot_ids=ids[id_idx],height=4,width=3)
})
layout <- paste0('\n',paste(layout_lines,collapse='\n'),'\n')

patches <- c(rows_titles,plot_list)
names(patches) <- ids

sbc_panel <- wrap_plots(patches, guides = 'collect', design = layout)
save_tikz_plot(sbc_panel, width = 6, filename = "tex/calibration_plot.tex")
