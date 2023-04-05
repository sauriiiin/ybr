##### YBR TUNI EXPERIMENT
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 03/27/2023

##### INITIALIZE
source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/initialize.R')
source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/fit_sum.R')

expt_info <- readxl::read_xlsx('/home/sbp29/RAW_Data/YBR_TUNI/Expt21/20230309_YBRs_INFO.xlsx')

##### GATHER DATA
fit_data <- NULL
for (e in unique(expt_info$expt_id)) {
  for (d in unique(expt_info$density[expt_info$expt_id == e])) {
    temp <- dbGetQuery(conn, sprintf('select * from %s_FS_ALL_%d_FITNESS a, %s_pos2coor b
                                     where a.pos = b.pos
                                     order by a.hours, b.plate_no, b.plate_col, b.plate_row', e, d, e))
    fit_data <- rbind(fit_data, temp[,-8])
  }
}
head(fit_data)
fit_data$rep <- as.numeric(str_trunc(as.character(fit_data$pos), 4, side = 'left', ellipsis = ''))
fit_data$condition[fit_data$density == 1536 & fit_data$plate_no %in% c(1,2)] <- 'YPDA (pH 6.8)'
fit_data$condition[fit_data$density == 1536 & fit_data$plate_no %in% c(3,4)] <- 'YPDA (pH 7.5)'
fit_data$condition[fit_data$density == 1536 & fit_data$plate_no %in% c(5,6)] <- 'DMSO (pH 7.5)'
fit_data$condition[fit_data$density == 1536 & fit_data$plate_no %in% c(7,8)] <- 'Tunicamycin 1ug/ml (pH 7.5)'

fit_data$condition[fit_data$density == 6144 & fit_data$plate_no %in% c(1)] <- 'YPDA (pH 6.8)'
fit_data$condition[fit_data$density == 6144 & fit_data$plate_no %in% c(2)] <- 'YPDA (pH 7.5)'
fit_data$condition[fit_data$density == 6144 & fit_data$plate_no %in% c(3)] <- 'DMSO (pH 7.5)'
fit_data$condition[fit_data$density == 6144 & fit_data$plate_no %in% c(4)] <- 'Tunicamycin 1ug/ml (pH 7.5)'

fit_data$condition <- factor(fit_data$condition, levels = c('YPDA (pH 6.8)', 'YPDA (pH 7.5)', 'DMSO (pH 7.5)', 'Tunicamycin 1ug/ml (pH 7.5)'))

fit_data %>%
  group_by(density, plate_no) %>%
  count()

##### CLEAN DATA AND GENERATE SUMMARY
fit_sum <- fit_data %>%
  group_by(density, condition, hours, rep, strain_id, orf_name) %>%
  summarise(cs.median = median(average, na.rm = T), cs.mad = mad(average, na.rm = T),
            fit.median = median(fitness, na.rm = T), fit.mad = mad(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

fit_data <- merge(fit_data, fit_sum,
                      by = c('density','condition','hours','rep','strain_id','orf_name'), all = T)

fit_data$average[fit_data$average < (fit_data$cs.median - 2*fit_data$cs.mad) |
                       fit_data$average > (fit_data$cs.median + 2*fit_data$cs.mad)] <- NA
fit_data$fitness[fit_data$fitness < (fit_data$fit.median - 2*fit_data$fit.mad) |
                   fit_data$fitness > (fit_data$fit.median + 2*fit_data$fit.mad)] <- NA

fit_sum <- fit_data %>%
  group_by(density, condition, hours, rep, strain_id, orf_name) %>%
  summarise(cs.median = median(average, na.rm = T),
            fit.median = median(fitness, na.rm = T),
            .groups = 'keep') %>%
  data.frame()



fit_data %>%
  filter(orf_name != 'NULL') %>%
  ggplot() +
  stat_summary(aes(x = hours, y = average, col = orf_name),
               fun=mean, geom="line", lwd = 1) +
  stat_summary(aes(x = hours, y = average, fill = orf_name),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  facet_grid(density ~ condition, scales = 'free_y') +
  labs(x="Hours", y="Colony Size (pix.)")+
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, margin = margin(0,0,0,0)),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = ggtext::element_markdown(size = txt),
        legend.position = 'bottom',
        legend.spacing = unit(0.1,"mm"),
        legend.key.size = unit(5, "mm"),
        legend.background = element_rect(fill = 'transparent'),
        strip.text = element_text(size = txt, margin = margin(0.1,0,0.1,0, "mm")))


fit_data %>%
  filter(orf_name != 'NULL', density == 6144, hours <= 24) %>%
  ggplot() +
  stat_summary(aes(x = hours, y = average, col = orf_name),
               fun=mean, geom="line", lwd = 1) +
  stat_summary(aes(x = hours, y = average, fill = orf_name),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  scale_y_log10() +
  facet_grid(density ~ condition, scales = 'free_y') +
  labs(x="Hours", y="Colony Size (pix.)")+
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, margin = margin(0,0,0,0)),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = ggtext::element_markdown(size = txt),
        legend.position = 'bottom',
        legend.spacing = unit(0.1,"mm"),
        legend.key.size = unit(5, "mm"),
        legend.background = element_rect(fill = 'transparent'),
        strip.text = element_text(size = txt, margin = margin(0.1,0,0.1,0, "mm")))

fit_data %>%
  filter(orf_name %notin% c('NULL','REF_BY4741'), density == 6144, hours == 24) %>%
  ggplot(aes(x = orf_name, y = average)) +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() +
  facet_grid(.~condition) 


fit_sum %>%
  filter(orf_name %notin% c('NULL','REF_BY4741'), density == 6144, hours == 24) %>%
  ggplot(aes(x = orf_name, y = fit.median)) +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() +
  facet_grid(.~condition)


merge(fit_sum %>% filter(density == 6144, hours == 24, condition == 'DMSO (pH 7.5)'),
      fit_sum %>% filter(density == 6144, hours == 24, condition == 'Tunicamycin 1ug/ml (pH 7.5)'),
      by = c('density','hours','rep','strain_id','orf_name'),
      suffixes = c('_DMSO','_TUNI')) %>%
  filter(orf_name %notin% c('NULL','REF_BY4741')) %>%
  ggplot(aes(x = fit.median_DMSO, y = fit.median_TUNI)) +
  geom_abline() +
  geom_point(aes(col = orf_name))

