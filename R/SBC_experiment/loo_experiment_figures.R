library(ggplot2)
library(readr)

mean_r2d2 <- 0.3
prec_r2d2 <- 20
conc_r2d2 <- 0.2

elpd_dat <- read_tsv(paste0('data/loo_experiment_r2d2_',mean_r2d2,'_',prec_r2d2,'_',conc_r2d2,'.tsv'))

p3 <- ggplot(elpd_dat) +
      geom_col(aes(as.factor(iter),-elpd_diff_normal_r2d2normal),fill='red',alpha=0.2,col='black') +
      geom_errorbar(aes(x=as.factor(iter),
                        ymin=-elpd_diff_normal_r2d2normal-1.96*se_diff_normal_r2d2normal,
                        ymax=-elpd_diff_normal_r2d2normal+1.96*se_diff_normal_r2d2normal),
                    width=0.2) +
      scale_y_continuous(expand=c(0,0)) +
      ylab(expression(Delta~elpd))+
      xlab('Iteration') +
      #ggtitle(paste0('mean=',mean_r2d2,', prec=',prec_r2d2,', conc=',conc_r2d2)) +
      theme_classic()

#Average elpd across iterations
-mean(elpd_dat$elpd_diff_normal_r2d2normal)
#Avergae standard error across iterations
mean(elpd_dat$se_diff_normal_r2d2normal)

p1 + p2 + p3
