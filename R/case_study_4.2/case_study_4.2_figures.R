library(brms)
library(projpred)
library(dplyr)
library(ggplot2)
library(gamm4)
library(scam)
library(geomtextpath)
library(tidyr)
source('R/aux/aux_plotting.R')

load('R/dat/case_study_4.2_R2_0.5_rho_0.1_fix_betaFALSE.RData')
sel_df <- vs_forward_loo_validated$summary
sel_df_not_null_mod <- filter(sel_df,size!=0)
mono.spline <- scam(diff/diff.se ~ s(size, k = 10, bs = "mpi", m = 2), data = sel_df_not_null_mod)

spline_df <- tibble(size=sel_df_not_null_mod$size[sel_df_not_null_mod$size>0],
                    mono.spline.diff = mono.spline$fit*sel_df_not_null_mod[,'diff.se'],
                    mono.spline.diff.se = sqrt(mono.spline$sig2)*sel_df_not_null_mod[,'diff.se'])
mod_size_candidates <- c('confidence_heuristic' = spline_df$size[min(which(spline_df$mono.spline.diff + spline_df$mono.spline.diff.se>=0))],
                         'elpd_difference_heuristic'=spline_df$size[min(which(spline_df$mono.spline.diff >= -4))])
elpd_ref <- slice(vs_forward_loo_validated$summary,1) %>% mutate(elpd_ref=elpd.loo-diff) %>% .$elpd_ref
legend_names <- c('Single-path forward selection','Cross-validated forward selection')
sel_df <- bind_rows(vs_forward_loo$summary,sel_df,.id='varsel') %>%
          mutate(varsel=factor(legend_names[as.integer(varsel)],levels=legend_names))
#### ELPD plot
elpd_plot <- ggplot() +
            geom_hline(yintercept=0,colour = "red", linetype = "longdash") +
            geom_line(data=filter(sel_df,varsel=='Single-path forward selection'),
                      aes(x=size,
                          y=diff),col='grey') +
            geom_pointrange(data=sel_df,
              aes(x=size,
                  y=diff,
                  ymin=diff-diff.se,
                  ymax=diff+diff.se,
                  col=varsel),
              position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0)
            ) +
            geom_ribbon(
              data = spline_df,
              aes(
                x = size,
                ymin = mono.spline.diff - mono.spline.diff.se,
                ymax = mono.spline.diff + mono.spline.diff.se
              ),
              alpha = 0.2, linetype = 2, fill = "blue"
            ) +
            geom_line(
              data = spline_df,
              aes(x = size,
                  y = mono.spline.diff
              ),
              col='blue') +
            geom_vline(xintercept = mod_size_candidates['confidence_heuristic'],
                       colour = "grey",
                       linetype = "longdash") +
            annotate("text",
                     x = 3,
                     y = 3,
                     colour = "red",
                     label = "Reference model elpd",size=3) +
            scale_x_continuous(breaks = seq(0,50,by=2)) +
            scale_color_manual(values=c('grey','black'),name='') +
            xlab('Model size') +
            ylab('Delta elpd') +
            theme_bw() +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  legend.position = c(0.8,0.45))

save_tikz_plot(elpd_plot,width=6.5,height=4,filename = 'tex/elpd_plot_case_study_4.2.tex')
