library(brms)
library(projpred)
library(dplyr)
library(ggplot2)
library(gamm4)
library(scam)
library(geomtextpath)
library(tidyr)
source('R/aux/aux_plotting.R')

load('results/varsel_objects_R2_0.5_rho_0.1_fix_betaFALSE.RData')
sel_df <- vs_forward_loo_validated$summary
mono.spline <- scam(diff/diff.se ~ s(size, k = 10, bs = "mpi", m = 2), data = sel_df)
spline_df <- tibble(size=sel_df$size,
                    mono.spline.diff = mono.spline$fit*sel_df[,'diff.se'],
                    mono.spline.diff.se = sqrt(mono.spline$sig2)*sel_df[,'diff.se'])
mod_size_candidates <- c('confidence_heuristic' = spline_df$size[min(which(spline_df$mono.spline.diff + spline_df$mono.spline.diff.se>=0))],
                         'elpd_difference_heuristic'=spline_df$size[min(which(spline_df$mono.spline.diff >= -4))])
elpd_ref <- slice(vs_forward_loo_validated$summary,1) %>% mutate(elpd_ref=elpd.loo-diff) %>% .$elpd_ref
sel_df <- bind_rows(list('Cross-validated forward selection'=sel_df,'Single-path forward selection'=vs_forward_loo$summary),.id='varsel')
#### ELPD plot
elpd_plot <- ggplot(data=sel_df) +
            geom_hline(yintercept=elpd_ref,colour = "red", linetype = "longdash") +
            geom_line(aes(x=size,
                          y=elpd.loo,
                          col=varsel,
                          alpha=varsel)) +
            geom_pointrange(
              aes(x=size,
                  y=elpd.loo,
                  ymin=elpd.loo-se,
                  ymax=elpd.loo+se,
                  col=varsel,
                  alpha=varsel),
              position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0)
            ) +
            geom_ribbon(
              data = spline_df,
              aes(
                x = size,
                ymin = elpd_ref+mono.spline.diff - mono.spline.diff.se,
                ymax = elpd_ref+mono.spline.diff + mono.spline.diff.se
              ),
              alpha = 0.2, linetype = 2, fill = "blue"
            ) +
            geom_line(
              data = spline_df,
              aes(x = size,
                  y = elpd_ref+mono.spline.diff
              ),
              col='blue') +
            geom_vline(xintercept = mod_size_candidates['confidence_heuristic'],
                       colour = "grey",
                       linetype = "longdash") +
            annotate("text",
                     x = 3,
                     y = elpd_ref + 3,
                     colour = "red",
                     label = "Reference model elpd",size=3) +
            scale_x_continuous(breaks = seq(0,50,by=2)) +
            scale_color_manual(values=c('blue','grey'),
                               name='') +
            scale_alpha_discrete(range = c(1, 0.3),
                                 name='') +
            xlab('Model size') +
            ylab('elpd') +
            theme_bw() +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  legend.position = c(0.8,0.55))

save_tikz_plot(elpd_plot,width=6.5,height=4,filename = 'tex/elpd_plot_case_study_4.2.tex')

######## Projected posteriors ######
extract_post <- function(mod,mod_name){
  post_mat <- as.matrix(mod)[,c('b_X2','b_X7','sigma')]
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
mod_types <- c('ref'='Reference model','prj'='Projected model')
var_names <- c('b_X2'='$beta_2$','b_X7'='$beta_7$','sigma'='$sigma$')
prj_submodel <- project(vs_forward_loo_validated$refmodel,solution_terms=vs_forward_loo_validated$solution_terms[seq_len(mod_size_candidates['confidence_heuristic'])],ndraws=4000)
mod_list <- list(vs_forward_loo_validated$refmodel$fit,prj_submodel)
posterior_dat <- bind_cols(extract_post(vs_forward_loo_validated$refmodel$fit,'ref_all'),
                           extract_post(prj_submodel,'prj_submodel')) %>%
  pivot_longer(cols=everything(),names_to='var',values_to='value') %>%
  separate(var,sep = ':',into = c('model','var'),remove = F) %>%
  separate(model,sep = '_',into = c('model_type','size'),remove = F) %>%
  mutate(mod_type_name=mod_types[model_type],
         var_name=var_names[var])

pp_dat <- bind_rows(extract_post_pred(vs_forward_loo_validated$refmodel$fit,mod_types['ref']),
                    extract_post_pred(prj_submodel,mod_types['prj'])) %>%
          mutate(var_name='$y$')
pp_summary <- bind_rows(extract_post_pred(vs_forward_loo_validated$refmodel$fit,mod_types['ref']),
                        extract_post_pred(prj_submodel,mod_types['prj'])) %>%
  group_by(mod_type_name,datapoint) %>%
  summarise(lower=quantile(value,0.025),
            mean = mean(value),
            upper=quantile(value,0.975)) %>%
  mutate(var_name='$y$')

cols <- c('#BC3C29FF','#0072B5FF','#E18727FF','#20854EFF','#7876B1FF','#6F99ADFF','#EE4C97FF')[seq_len(2)]
plot_dat <- select(pp_summary,mod_type_name,var_name,value=mean)
prj_post_plot <- ggplot(plot_dat) +
  geom_density(data=mutate(vs_forward_loo_validated$refmodel$fetch_data(),var_name='$y$'),aes(y),color=NA,fill='black',alpha=0.2,adjust=1) +
  stat_density(aes(value,col=mod_type_name),geom='line',position = "identity",adjust=1,size=0.8) +
  facet_wrap(~var_name,nrow=1,scale='free') +
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

save_tikz_plot(prj_post_plot,width=6,height=3,filename = 'tex/projected_posteriors_case_study_4.2.tex')
