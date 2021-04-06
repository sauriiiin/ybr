library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)
library(egg)
library(zoo)
library(ggrepel)
library(plotly)
library(scales)
library(reshape2)
library(RMariaDB)
library(RColorBrewer)

source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

out_path <- 'output/figs';

remove_outliers <- function(x) {
  qnt <- quantile(x[!is.na(x)], probs=c(.05, .95), na.rm = TRUE)
  y <- x
  y[x < qnt[1]] <- NA
  y[x > qnt[2]] <- NA
  y
}

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9

##### PLATEMAPS
ybrs_pm <- dbGetQuery(conn, 'select * from YBRs_pos2coor a, YBRs_pos2orf_name b
                  where a.pos = b.pos')

ggplot(ybrs_pm[ybrs_pm$density == 1536,],
       aes(x = col, y = row, fill = orf_name)) +
  geom_tile(col = 'black') +
  scale_y_reverse() +
  facet_wrap(.~density*plate, 
             scale = 'free',
             ncol = 2) +
  theme_linedraw()

plyr::count(ybrs_pm, vars = c('density','orf_name'))

##### FITNESS RESULTS
fitness <- dbGetQuery(conn, 'select *
                   from YBRs_YPD_1536_FITNESS
                   union
                   select *
                   from YBRs_YPD_6144_FITNESS')
fitness$density[fitness$pos <= 999999] <- 1536
fitness$density[fitness$pos >= 999999] <- 6144
plyr::count(fitness, vars = 'density')

fitness$orf_name <- factor(fitness$orf_name, levels = c('REF',
                                                        'FY4',
                                                        'BY4741',
                                                        'BY4741_KAN',
                                                        'BY4742',
                                                        'BY4742_KAN',
                                                        'CRISPY',
                                                        'ATG',
                                                        'PAM',
                                                        'TATA',
                                                        'C_STOP',
                                                        'NC_STOP',
                                                        'NC_STOP2',
                                                        'BOR'))


# fit.stats <- dbGetQuery(conn, 'select a.orf_name, a.cs_mean 1536_mean, a.cs_median 1536_median, a.cs_std 1536_std,
#                         b.cs_mean 6144_mean, b.cs_median 6144_median, b.cs_std 6144_std
#                         from YBRs_YPD_1536_FITNESS_STATS a, YBRs_YPD_6144_FITNESS_STATS b
#                         where a.orf_name = b.orf_name')
# 
# fit.stats$orf_name <- factor(fit.stats$orf_name, levels = c('REF',
#                                                         'FY4',
#                                                         'BY4741',
#                                                         'BY4741_KAN',
#                                                         'BY4742',
#                                                         'BY4742_KAN',
#                                                         'CRISPY',
#                                                         'ATG',
#                                                         'PAM',
#                                                         'TATA',
#                                                         'C_STOP',
#                                                         'NC_STOP',
#                                                         'NC_STOP2',
#                                                         'BOR'))

## CLEAN fitness
for (d in unique(fitness$density)) {
  for (h in unique(fitness$hours[fitness$density == d])) {
    for (o in unique(fitness$orf_name[fitness$density == d & fitness$hours == h])) {
      fitness$fitness[fitness$density == d & fitness$hours == h & fitness$orf_name == o] <-
        remove_outliers(fitness$fitness[fitness$density == d & fitness$hours == h & fitness$orf_name == o])
    }
  }
}

## PLOT fitness
ggplot(fitness[fitness$orf_name != 'BOR' & fitness$hours > 5,],
       aes(x = orf_name, y = fitness, fill = orf_name)) +
  geom_hline(yintercept = 1, col = 'red') +
  geom_jitter(size = 0.05, alpha = 0.7) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5) +
  labs(title = 'YPD') +
  scale_fill_discrete(guide = F) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5)) +
  facet_wrap(.~density*hours) +
  theme(plot.title = element_text(size = titles,
                                  face = 'bold',
                                  hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        strip.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"))

ggsave(sprintf("%sFITNESS_YPD.jpg",out_path),
       height = 150, width = two.c, units = 'mm',
       dpi = 300)

# ggplot(fit.stats[fit.stats$orf_name != 'BOR',],
#        aes(x = `1536_median`, y = `6144_median`, col = orf_name)) +
#   geom_point() +
#   geom_text(aes(label = orf_name))

##### YBRs3 RESULTS
##### PLATEMAPS
ybrs3_pm <- dbGetQuery(conn, 'select * from YBRs3_pos2coor a, YBRs3_pos2orf_name b
                  where a.pos = b.pos')
ybrs3_pm$orf_name[ybrs3_pm$orf_name == 'STOP_NC2'] <- 'STOP_NC'
ybrs3_pm$orf_name <- factor(ybrs3_pm$orf_name, levels = c('REF','PAM','crispy','ATG','TATA','STOP_C','STOP_NC','ire1','hog1',
                                                          'npp2','ptp2','uth1','mid2','stl1','ssm4','trp1','trp3','FY4','BY4741',
                                                          'BY4742','kan41','kan42','aft1','fet12','BOR'))

colourCount = length(unique(ybrs3_pm$orf_name))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

for (d in unique(ybrs3_pm$density)) {
  for (p in unique(ybrs3_pm$plate[ybrs3_pm$density == d])) {
    nam <- sprintf('d%dp%d',d,p)
    
    temp <- ggplot(ybrs3_pm[ybrs3_pm$density == d & ybrs3_pm$plate == p,],
                   aes(x = col, y = row, fill = orf_name)) +
      geom_tile(col = 'black') +
      labs(title = sprintf('%d Density - Plate %d', d, p),
           x = 'Column',
           y = 'Row') +
      scale_y_reverse() +
      scale_fill_manual(name = 'Strain',
                        breaks = c('REF','PAM','crispy','ATG','TATA','STOP_C','STOP_NC','ire1','hog1',
                                   'npp2','ptp2','uth1','mid2','stl1','ssm4','trp1','trp3','FY4','BY4741',
                                   'BY4742','kan41','kan42','aft1','fet12','BOR'),
                        values = getPalette(colourCount),
                        drop = F) +
      # facet_wrap(.~density*plate, 
      #            scale = 'free',
      #            ncol = 2) +
      theme_linedraw() +
      theme(plot.title = element_text(size = titles,
                                      face = 'bold',
                                      hjust = 0.5),
            axis.title = element_text(size = titles),
            axis.text = element_text(size = txt),
            legend.title = element_text(size = titles),
            legend.text = element_text(size = txt),
            legend.position = 'bottom',
            legend.key.size = unit(3, "mm"),
            legend.box.spacing = unit(0.5,"mm"),
            strip.text = element_text(size = txt,
                                      face = 'bold',
                                      margin = margin(0.1,0,0.1,0, "mm"))) +
      guides(fill = guide_legend(nrow=4))
    
    assign(nam, temp)
  }
}

ggpubr::ggarrange(d96p1, d96p2, d96p3,
          d384p1, d384p2, NULL,
          d1536p1, d1536p2, NULL,
          common.legend = T,legend = 'bottom')
ggsave(sprintf("%sPM_YBRs3.jpg",out_path),
       height = two.c*1.5, width = two.c*2, units = 'mm',
       dpi = 300)

#### FITNESS FS1
expt_names <- data.frame(id = c('CTRL','A','B','C','D','E','F','G1','G2','H','I','J1','J2','K','L1','L2','M','N'),
                         condition = c('YPD','YPD @ 37C','YPD @ 20C','SC + GLU','YPD + 6% EtOH','YPG (3% Glyerol)','YPD + 1M Sorbitol',
  'YPD + 0.05% CW','YPD + 0.1% CW','SC + GLU @ 50C','YPD + 10mM Caffeine','YPD + 1M NaCl','YPD + 1.5M NaCl',
  'YPD + 0.1% SDS','YPD + 1ug/ml Tunicamycin','YPD + DMSO','YPD + 1ug/ml Congo Red','SC + CASE + GAL'))

ybrs3_fit <- NULL
for (e in expt_names$id) {
  temp <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate
                                   from YBRs3_%s_1536_FITNESS a, YBRs3_pos2coor b
                                   where a.pos = b.pos and a.hours > 16', e))
  temp$id <- e
  temp$condition <- expt_names$condition[expt_names$id == e]
  for (h in unique(temp$hours)) {
    for (o in unique(temp$orf_name[temp$hours == h])) {
      temp$fitness[temp$hours == h & temp$orf_name == o] <-
        remove_outliers(temp$fitness[temp$hours == h & temp$orf_name == o])
    }
  }
  ybrs3_fit <- rbind(ybrs3_fit, temp)
}

ybrs3_fit$orf_name[ybrs3_fit$orf_name == 'STOP_NC2'] <- 'STOP_NC'

ybrs3_fit$orf_name <- factor(ybrs3_fit$orf_name, levels = c('REF','PAM','crispy','ATG','TATA','STOP_C','STOP_NC','ire1','hog1',
                                                            'npp2','ptp2','uth1','mid2','stl1','ssm4','trp1','trp3','FY4','BY4741',
                                                            'BY4742','kan41','kan42','aft1','fet12','BOR'))

ybrs3_fit$condition <- factor(ybrs3_fit$condition , levels = c('YPD','YPD @ 37C','YPD @ 20C','SC + GLU','YPD + 6% EtOH','YPG (3% Glyerol)','YPD + 1M Sorbitol',
                                                               'YPD + 0.05% CW','YPD + 0.1% CW','SC + GLU @ 50C','YPD + 10mM Caffeine','YPD + 1M NaCl','YPD + 1.5M NaCl',
                                                               'YPD + 0.1% SDS','YPD + 1ug/ml Tunicamycin','YPD + DMSO','YPD + 1ug/ml Congo Red','SC + CASE + GAL'))

ybrs3_fit$fitness[ybrs3_fit$id == 'J1' & ybrs3_fit$plate == 1] <- NA
ybrs3_fit$fitness[ybrs3_fit$id == 'J2' & ybrs3_fit$plate == 2] <- NA
ybrs3_fit$fitness[ybrs3_fit$id == 'E' & ybrs3_fit$plate == 1] <- NA
ybrs3_fit <- ybrs3_fit[ybrs3_fit$orf_name != 'BOR',]


# ybrs3.fit.plot <- ggplot(ybrs3_fit[!is.na(ybrs3_fit$condition),],
#        aes(x = orf_name, y = fitness, fill = orf_name)) +
#   geom_hline(yintercept = 1, col = 'red') +
#   # geom_jitter(size = 0.05, alpha = 0.7) +
#   # geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5) +
#   geom_boxplot(outlier.shape = NA) +
#   stat_compare_means(method = 'wilcox.test', ref.group = "REF",
#                      label = "p.signif", label.y = 1.7,
#                      size = 1.5, angle = 90) +
#   scale_fill_discrete(guide = F) +
#   facet_wrap(.~condition, nrow = 3) +
#   coord_cartesian(ylim = c(0.25,1.75)) +
#   theme_linedraw() +
#   theme(plot.title = element_text(size = titles,
#                                   face = 'bold',
#                                   hjust = 0.5),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         axis.text.x = element_text(angle = 90,
#                                    vjust = 0.5),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.position = 'bottom',
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   face = 'bold',
#                                   margin = margin(0.1,0,0.1,0, "mm")))
# 
# ggsave(sprintf("%sFITNESS_ALL_CLEAN.jpg",out_path),ybrs3.fit.plot,
#        height = 150, width = two.c*2, units = 'mm',
#        dpi = 300)

#### FITNESS FS1 REDO
ybrs3_r_fit <- NULL
for (e in c('CTRL','J1','J2')) {
  temp <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate
                                   from YBRs3_R_%s_1536_FITNESS a, YBRs3_pos2coor b
                                   where a.pos = b.pos and a.hours > 16', e))
  temp$id <- e
  temp$condition <- expt_names$condition[expt_names$id == e]
  for (h in unique(temp$hours)) {
    for (o in unique(temp$orf_name[temp$hours == h])) {
      temp$fitness[temp$hours == h & temp$orf_name == o] <-
        remove_outliers(temp$fitness[temp$hours == h & temp$orf_name == o])
    }
  }
  ybrs3_r_fit <- rbind(ybrs3_r_fit, temp)
}

ybrs3_r_fit$orf_name[ybrs3_r_fit$orf_name == 'STOP_NC2'] <- 'STOP_NC'
ybrs3_r_fit$orf_name <- factor(ybrs3_r_fit$orf_name, levels = c('REF','PAM','crispy','ATG','TATA','STOP_C','STOP_NC','ire1','hog1',
                                                            'npp2','ptp2','uth1','mid2','stl1','ssm4','trp1','trp3','FY4','BY4741',
                                                            'BY4742','kan41','kan42','aft1','fet12','BOR'))

ybrs3_r_fit$condition <- factor(ybrs3_r_fit$condition , levels = c('YPD','YPD @ 37C','YPD @ 20C','SC + GLU','YPD + 6% EtOH','YPG (3% Glyerol)','YPD + 1M Sorbitol',
                                                               'YPD + 0.05% CW','YPD + 0.1% CW','SC + GLU @ 50C','YPD + 10mM Caffeine','YPD + 1M NaCl','YPD + 1.5M NaCl',
                                                               'YPD + 0.1% SDS','YPD + 1ug/ml Tunicamycin','YPD + DMSO','YPD + 1ug/ml Congo Red','SC + CASE + GAL'))

ybrs3_r_fit <- ybrs3_r_fit[ybrs3_r_fit$orf_name != 'BOR',]

ybrs3_r_fit$attempt <- 2
ybrs3_fit$attempt <- 1

ybrs3_all_fit <- rbind(ybrs3_fit, ybrs3_r_fit)

ybrs3.r.fit.plot <- ggplot(ybrs3_r_fit[!is.na(ybrs3_r_fit$condition),],
       aes(x = orf_name, y = fitness, fill = orf_name)) +
  geom_hline(yintercept = 1, col = 'red') +
  # geom_jitter(size = 0.05, alpha = 0.7) +
  # geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = 'wilcox.test', ref.group = "REF",
                     label = "p.signif", label.y = 1.7,
                     size = 1.5, angle = 90) +
  scale_fill_discrete(guide = F) +
  facet_wrap(.~condition, nrow = 3) +
  coord_cartesian(ylim = c(0.25,1.75)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles,
                                  face = 'bold',
                                  hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))

ybrs3.all.fit.plot <- ggplot(ybrs3_all_fit[!is.na(ybrs3_all_fit$condition),],
                           aes(x = orf_name, y = fitness, fill = orf_name)) +
  geom_hline(yintercept = 1, col = 'red') +
  # geom_jitter(size = 0.05, alpha = 0.7) +
  # geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = 'wilcox.test', ref.group = "REF",
                     label = "p.signif", label.y = 1.7,
                     size = 1.5, angle = 90) +
  scale_fill_discrete(guide = F) +
  facet_wrap(.~attempt*condition, nrow = 3) +
  coord_cartesian(ylim = c(0.25,1.75)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles,
                                  face = 'bold',
                                  hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))

ggsave(sprintf("%sFITNESS_ALL_CLEAN_wREDO.jpg",out_path),
       height = 150, width = two.c*2 + two.c*2/6, units = 'mm',
       dpi = 300)

##### REF AND YBR MUTANTS ONLY

ybrs3.fit.plot <- ggplot(ybrs3_all_fit[!is.na(ybrs3_all_fit$condition) &
                                         ybrs3_all_fit$orf_name %in% c('REF','PAM','crispy','ATG','TATA','STOP_C','STOP_NC'),],
                             aes(x = orf_name, y = fitness, fill = orf_name)) +
  geom_hline(yintercept = 1, col = 'red') +
  # geom_jitter(size = 0.05, alpha = 0.7) +
  # geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = 'wilcox.test', ref.group = "REF",
                     label = "p.signif", label.y = 1.15,
                     size = 1.5, angle = 90, vjust = 0.5) +
  scale_fill_discrete(guide = F) +
  facet_wrap(.~attempt*condition, nrow = 3) +
  coord_cartesian(ylim = c(0.85,1.15)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles,
                                  face = 'bold',
                                  hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))

ggsave(sprintf("%sFITNESS_YBRs_CLEAN_wREDO.jpg",out_path), ybrs3.fit.plot,
       height = 150, width = two.c*2 + two.c*2/6, units = 'mm',
       dpi = 300)


##### CONDITION WISE with CONTROLS
main.mutants <- c('REF','PAM','crispy','ATG','TATA','STOP_C','STOP_NC')

id <- 'M'
ggplot(ybrs3_all_fit[ybrs3_all_fit$condition == expt_names$condition[expt_names$id == id] &
                       ybrs3_all_fit$orf_name %in% c(main.mutants, c('mid2','stl1')),],
       aes(x = orf_name, y = fitness, fill = orf_name)) +
  geom_hline(yintercept = 1, col = 'red') +
  # geom_jitter(size = 0.05, alpha = 0.7) +
  # geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = 'wilcox.test', ref.group = "REF",
                     label = "p.signif", label.y = 1.3,
                     size = 1.5, angle = 90, vjust = 0.5) +
  scale_fill_discrete(guide = F) +
  labs(title = expt_names$condition[expt_names$id == id],
       x = 'Mutant',
       y = 'Fitness') +
  # facet_wrap(nrow = 3) +
  coord_cartesian(ylim = c(0.7,1.3)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles,
                                  face = 'bold',
                                  hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))

ggsave(sprintf("%sFITNESS_YBRs_%s.jpg",out_path,id),
       height = one.c, width = one.c, units = 'mm',
       dpi = 300)
