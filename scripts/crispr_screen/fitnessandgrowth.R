
source('/home/sbp29/R/Projects/ybr/scripts/crispr_screen/initialize.R')

expt.file <- read_xlsx(path = '/home/sbp29/RAW_Data/YBR_CRISPR/YBR_CRISPR_INFO.xlsx') %>% data.frame()
head(expt.file)

failed.pcr <- c(20000472, 20000473, 20000482, 20000486, 20000497, 20000498, 20000499, 20000500,
                20000502, 20000503, 20000508, 20000538, 20000539, 20000541, 20000547, 20000548,
                20000549, 20000554)

conditions <- c('YPDA','SA','HU','HO',
                'DM','TN','FL')
mutants <- c('REF','PAM1','PAM1_ATG12','PAM1_STOP1C','PAM1_STOP13C',
             'WT','STOP23','STOP23C','STOP1C','STOP2C','STOP13C','DEL','KANMX',
             'BOR','NULL')

ybr_crispr_fitness <- NULL
for (e in unique(expt.file$expt_id)) {
  # for (s in unique(expt.file$stage_id[expt.file$expt_id == e])) {
  s <- 'FS'
    for (a in unique(expt.file$arm[expt.file$expt_id == e & expt.file$stage_id == s])) {
      for (d in unique(expt.file$density[expt.file$expt_id == e & expt.file$stage_id == s & expt.file$arm == a])) {
        temp.fit <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate_no, b.plate_col, b.plate_row
                                         from %s_%s_%s_%d_FITNESS a, %s_pos2coor b
                                         where a.pos = b.pos
                                         order by a.hours, b.plate_no, b.plate_col, b.plate_row', e, s, a, d, e))
        temp.fit$expt_id <- e
        temp.fit$stage_id <- s
        temp.fit$condition <- a
        
        ybr_crispr_fitness <- rbind(ybr_crispr_fitness, temp.fit)
      }
    }
  # }
}
ybr_crispr_fitness$fitness[ybr_crispr_fitness$density == 384 &
                             ((ybr_crispr_fitness$plate_col %in% c(22:24) &
                             ybr_crispr_fitness$plate_row %in% c(8,15)) |
                               (ybr_crispr_fitness$plate_col == 22 &
                                  ybr_crispr_fitness$plate_row %in% c(7:14)))] <- NA
ybr_crispr_fitness$average[ybr_crispr_fitness$density == 384 &
                             ((ybr_crispr_fitness$plate_col %in% c(22:24) &
                                 ybr_crispr_fitness$plate_row %in% c(8,15)) |
                                (ybr_crispr_fitness$plate_col == 22 &
                                   ybr_crispr_fitness$plate_row %in% c(7:14)))] <- NA

ybr_crispr_fitness$fitness[ybr_crispr_fitness$density == 1536 &
                             ((ybr_crispr_fitness$plate_col %in% c(44:48) &
                                 ybr_crispr_fitness$plate_row %in% c(16,29)) |
                                (ybr_crispr_fitness$plate_col == 44 &
                                   ybr_crispr_fitness$plate_row %in% c(17:28)))] <- NA
ybr_crispr_fitness$average[ybr_crispr_fitness$density == 1536 &
                             ((ybr_crispr_fitness$plate_col %in% c(44:48) &
                                 ybr_crispr_fitness$plate_row %in% c(16,29)) |
                                (ybr_crispr_fitness$plate_col == 44 &
                                   ybr_crispr_fitness$plate_row %in% c(17:28)))] <- NA

ybr_crispr_fitness$fitness[ybr_crispr_fitness$density == 6144 &
                             ((ybr_crispr_fitness$plate_col %in% c(88:96) &
                                 ybr_crispr_fitness$plate_row %in% c(32,57)) |
                                (ybr_crispr_fitness$plate_col == 88 &
                                   ybr_crispr_fitness$plate_row %in% c(33:56)))] <- NA
ybr_crispr_fitness$average[ybr_crispr_fitness$density == 6144 &
                             ((ybr_crispr_fitness$plate_col %in% c(88:96) &
                                 ybr_crispr_fitness$plate_row %in% c(32,57)) |
                                (ybr_crispr_fitness$plate_col == 88 &
                                   ybr_crispr_fitness$plate_row %in% c(33:56)))] <- NA



ybr_crispr_summary <- ybr_crispr_fitness %>%
  filter(orf_name %notin% c('BOR','NULL')) %>%
  group_by(expt_id, stage_id, density, condition, hours, strain_id, orf_name) %>%
  summarize(avg.median = median(average, na.rm = T), avg.mad = mad(average, na.rm = T),
            fitness.median = median(fitness, na.rm = T), fitness.mad = mad(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

##### REMOVE OUTLIERS
ybr_crispr_fitness <- merge(ybr_crispr_fitness, ybr_crispr_summary,
                      by = c('expt_id','stage_id','density','condition','hours','strain_id','orf_name'), all = T)

ybr_crispr_fitness$average[ybr_crispr_fitness$average < (ybr_crispr_fitness$avg.median - 2*ybr_crispr_fitness$avg.mad) |
                       ybr_crispr_fitness$average > (ybr_crispr_fitness$avg.median + 2*ybr_crispr_fitness$avg.mad)] <- NA
ybr_crispr_fitness$fitness[ybr_crispr_fitness$fitness < (ybr_crispr_fitness$fitness.median - 2*ybr_crispr_fitness$fitness.mad) |
                       ybr_crispr_fitness$fitness > (ybr_crispr_fitness$fitness.median + 2*ybr_crispr_fitness$fitness.mad)] <- NA
ybr_crispr_fitness <- ybr_crispr_fitness[,c(1:14)]

ybr_crispr_summary <- ybr_crispr_fitness %>%
  filter(orf_name %notin% c('BOR','NULL')) %>%
  group_by(expt_id, density, stage_id, condition, hours, strain_id, orf_name) %>%
  summarize(avg.median = median(average, na.rm = T), avg.mad = mad(average, na.rm = T),
            fitness.median = median(fitness, na.rm = T), fitness.mad = mad(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()
head(ybr_crispr_summary)

ybr_crispr_fitness$condition <- factor(ybr_crispr_fitness$condition, levels = conditions)
ybr_crispr_fitness$orf_name <- factor(ybr_crispr_fitness$orf_name, levels = mutants)
ybr_crispr_summary$condition <- factor(ybr_crispr_summary$condition, levels = conditions)
ybr_crispr_summary$orf_name <- factor(ybr_crispr_summary$orf_name, levels = mutants)

##### GROWTH CURVES
ybr_crispr_fitness %>%
  filter(strain_id %notin% failed.pcr,
         orf_name %notin% c('REF','BOR','NULL')) %>%
  ggplot(aes(x = hours, y = average)) +
  stat_summary(aes(fill = orf_name, group = orf_name),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.4) +
  stat_summary(aes(col = orf_name, group = orf_name),
               fun=mean, geom="line", lwd =1) +
  scale_color_discrete(guide = 'none') +
  scale_fill_discrete(guide = 'none') +
  facet_grid(density~condition,
             scales = 'free_y') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) 

ybr_crispr_summary %>%
  filter(strain_id %notin% failed.pcr,
         orf_name %notin% c('REF','BOR','NULL')) %>%
  ggplot(aes(x = hours, y = fitness.mad)) +
  stat_summary(aes(fill = orf_name, group = orf_name),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.4) +
  stat_summary(aes(col = orf_name, group = orf_name),
               fun=mean, geom="line", lwd =1) +
  scale_color_discrete(guide = 'none') +
  scale_fill_discrete(guide = 'none') +
  facet_grid(density~condition,
             scales = 'free_y') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) 

unique(ybr_crispr_fitness$hours)
ybr_crispr_fitness$saturation[ybr_crispr_fitness$density == 384] <- 42
ybr_crispr_fitness$saturation[ybr_crispr_fitness$density == 1536] <- 42 
ybr_crispr_fitness$saturation[ybr_crispr_fitness$density == 6144 &
                                ybr_crispr_fitness$condition %notin% c('FL','TN','SA')] <- 18
ybr_crispr_fitness$saturation[ybr_crispr_fitness$density == 6144 &
                                ybr_crispr_fitness$condition %in% c('FL','TN')] <- 24
ybr_crispr_fitness$saturation[ybr_crispr_fitness$density == 6144 &
                                ybr_crispr_fitness$condition %in% c('SA')] <- 42

ybr_crispr_summary$saturation[ybr_crispr_summary$density == 384] <- 42
ybr_crispr_summary$saturation[ybr_crispr_summary$density == 1536] <- 42 
ybr_crispr_summary$saturation[ybr_crispr_summary$density == 6144 &
                                ybr_crispr_summary$condition %notin% c('FL','TN','SA')] <- 18
ybr_crispr_summary$saturation[ybr_crispr_summary$density == 6144 &
                                ybr_crispr_summary$condition %in% c('FL','TN')] <- 24
ybr_crispr_summary$saturation[ybr_crispr_summary$density == 6144 &
                                ybr_crispr_summary$condition %in% c('SA')] <- 42

##### FITNESS PLOTS
plot.fit.sum <- ybr_crispr_fitness %>%
  filter(strain_id %notin% failed.pcr, hours == saturation,
         orf_name %notin% c('BOR','NULL','PAM1_STOP1C','PAM1_STOP13C',
                            'STOP1C','STOP13C')) %>%
  ggplot(aes(x = orf_name, y = fitness)) +
  # geom_jitter(aes(col = as.factor(strain_id)), size = 0.3) +
  geom_hline(yintercept = c(0.99,1.01), col = 'red', linetype = 'dashed', size = 0.2) +
  geom_hline(yintercept = c(0.95,1.05), col = 'red', linetype = 'solid', size = 0.2) +
  geom_boxplot(aes(fill = orf_name), outlier.shape = NA, size = 0.3) +
  scale_x_discrete(breaks = c('REF','PAM1','PAM1_ATG12',
                              'WT','STOP23','STOP23C','DEL','KANMX'),
                   labels = c('REF' = 'BY4741',
                              'PAM1' = 'PAM Mutant',
                              'PAM1_ATG12' = 'PAM + Two ATG Mutations',
                              'WT' = 'Wild-Type',
                              'STOP23' = 'Two Premature Stops',
                              'STOP23C' = 'Synonymous Controls for Stops',
                              'DEL' = 'CRISPR Deletion Mutant',
                              'KANMX' = 'KanMX Deletion Mutant')) +
  scale_color_discrete(guide = 'none') +
  scale_fill_discrete(guide = 'none') +
  labs(y = 'Fitness') +
  facet_grid(density~condition,
             labeller = labeller(condition = c('SA' = '1M NaCl',
                                                      'HU' = '100mM Hydroxyurea',
                                                      'HO' = '30mM Hydrogen Peroxide',
                                                      'FL' = '25ug/ml Fluconazole',
                                                      'TN' = '0.6uM Tunicamycin',
                                                      'YPDA' = 'YPDA',
                                                      'DM' = 'DMSO'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0.5,1.5))
ggsave(sprintf("%s/YBR_CRISPR_FITNESS_SUMMARY.jpg",fig_path), plot.fit.sum,
       height = two.c, width = two.c*2, units = 'mm',
       bg = 'white',
       dpi = 600)


ybr_crispr_fitness %>%
  filter(strain_id %notin% failed.pcr, hours == saturation,
         orf_name %in% c('PAM1','PAM1_ATG12',
                         'WT','STOP23','STOP23C','DEL',
                         'KANMX')) %>%
  ggplot(aes(y = orf_name, x = fitness)) +
  geom_density_ridges() +
  scale_color_discrete(guide = 'none') +
  facet_grid(density~condition) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.title.y = element_blank(),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) 


##### DIFFERENTIAL FITNESS
ybr_crispr_diff_fitness <- rbind(merge(ybr_crispr_fitness[ybr_crispr_fitness$condition %in% c('SA','HU','HO') &
                                                            ybr_crispr_fitness$hours == ybr_crispr_fitness$saturation,
                                                          c(3,4,6,7,8,10,11)], 
                                       ybr_crispr_fitness[ybr_crispr_fitness$condition %in% c('YPDA') &
                                                            ybr_crispr_fitness$hours == ybr_crispr_fitness$saturation,
                                                          c(3,4,6,7,8,10,11)],
                                       by = c('density','strain_id','orf_name','pos'),
                                       suffixes = c('_stress','_control')),
                                 merge(ybr_crispr_fitness[ybr_crispr_fitness$condition %in% c('TN','FL') &
                                                            ybr_crispr_fitness$hours == ybr_crispr_fitness$saturation,
                                                          c(3,4,6,7,8,10,11)], 
                                       ybr_crispr_fitness[ybr_crispr_fitness$condition %in% c('DM') &
                                                            ybr_crispr_fitness$hours == ybr_crispr_fitness$saturation,
                                                          c(3,4,6,7,8,10,11)],
                                       by = c('density','strain_id','orf_name','pos'),
                                       suffixes = c('_stress','_control')))

ybr_crispr_diff_fitness$fitness_diff <- ybr_crispr_diff_fitness$fitness_stress - ybr_crispr_diff_fitness$fitness_control


plot.diff.fit.sum <- ybr_crispr_diff_fitness %>%
  filter(strain_id %notin% failed.pcr,
         orf_name %notin% c('BOR','NULL','PAM1_STOP1C','PAM1_STOP13C',
                            'STOP1C','STOP13C')) %>%
  ggplot(aes(x = orf_name, y = fitness_diff)) +
  # geom_jitter(aes(col = as.factor(strain_id)), size = 0.3) +
  geom_hline(yintercept = c(-0.01,0.01), col = 'red', linetype = 'dashed', size = 0.2) +
  geom_hline(yintercept = c(-0.05,0.05), col = 'red', linetype = 'solid', size = 0.2) +
  geom_boxplot(aes(fill = orf_name), outlier.shape = NA, size = 0.3) +
  scale_color_discrete(guide = 'none') +
  scale_fill_discrete(guide = 'none') +
  scale_x_discrete(breaks = c('REF','PAM1','PAM1_ATG12',
                              'WT','STOP23','STOP23C','DEL','KANMX'),
                   labels = c('REF' = 'BY4741',
                              'PAM1' = 'PAM Mutant',
                              'PAM1_ATG12' = 'PAM + Two ATG Mutations',
                              'WT' = 'Wild-Type',
                              'STOP23' = 'Two Premature Stops',
                              'STOP23C' = 'Synonymous Controls for Stops',
                              'DEL' = 'CRISPR Deletion Mutant',
                              'KANMX' = 'KanMX Deletion Mutant')) +
  labs(y = 'Differential Fitness\n(Fitness in Stress - Control Condition)') +
  facet_grid(density~condition_stress,
             labeller = labeller(condition_stress = c('SA' = 'Salt',
                                                      'HU' = 'Hydroxyurea',
                                                      'HO' = 'Hydrogen Peroxide',
                                                      'FL' = 'Fluconazole',
                                                      'TN' = 'Tunicamycin'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) #+
  # coord_cartesian(ylim = c(-0.2,0.2))
ggsave(sprintf("%s/YBR_CRISPR_DIFF_FITNESS_SUMMARY.jpg",fig_path), plot.diff.fit.sum,
       height = two.c, width = two.c*2, units = 'mm',
       bg = 'white',
       dpi = 600)


ybr_crispr_diff_fitness %>%
  filter(strain_id %notin% failed.pcr,
         orf_name %notin% c('REF','BOR','NULL')) %>%
  ggplot(aes(x = fitness_stress, y = fitness_control)) +
  geom_point() +
  facet_grid(condition_stress~density) +
  coord_cartesian(xlim = c(0.5,1.5),
                  ylim = c(0.5,1.5))

