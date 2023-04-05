##### YBR SOLID SCREENS WITH NACL AND CACL2 (EXPT 45)
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 12/22/2022

##### INITIALIZE
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(dplyr)
library(stringr)
library(RMariaDB)

source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

`%notin%` <- Negate(`%in%`)

file_info <- readxl::read_xlsx(path = '/home/sbp29/RAW_Data/YBR_NACL_CACL/YBR_NACL_CACL_INFO.xlsx') %>% data.frame()
head(file_info)

arm_ids <- c("pH5", "NACL", "CACL2_60mM", "CACL2_60nM", "pH5_OLD")
cnd_ids <- c('YPD (pH 5.0)', '1M NaCl (pH 5.0)', '60mM CaCl2 (pH 7.5)', '60nM CaCl2 (pH 7.5)', 'YPD (pH 5.0) Old')
orf_ids <- c("BY4741_REF", "BY4741", "VMA3", "VPS16", "YBR_CRISP", "NULL")
strn_ids <- c("BY4741 (Ref.)", "BY4741","*vma3Δ::KanMX*","*vps16Δ::KanMX*","*ybr196c-aΔ*")

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

##### GATHER DATA
fit_data <- NULL
for (e in unique(file_info$expt_id)) {
  for (d in unique(file_info$density[file_info$expt_id == e])) {
    for (a in unique(file_info$arm[file_info$expt_id == e & file_info$density == d])) {
      for (s in unique(file_info$stage_id[file_info$expt_id == e & file_info$density == d & file_info$arm == a])) {
        temp <- dbGetQuery(conn, sprintf('select * from %s_%s_%s_%d_FITNESS', e, s, a, d))
        temp$expt_id <- e
        temp$density <- d
        temp$arm <- a
        temp$stage_id <- s
        
        fit_data <- rbind(fit_data, temp)
      }
    }
  }
}
fit_data$rep <- as.numeric(str_trunc(as.character(fit_data$pos), 4, side = 'left', ellipsis = ''))
fit_data$condition[fit_data$arm == 'pH5'] <- 'YPD (pH 5.0)'
fit_data$condition[fit_data$arm == 'pH5_OLD'] <- 'YPD (pH 5.0) Old'
fit_data$condition[fit_data$arm == 'NACL'] <- '1M NaCl (pH 5.0)'
fit_data$condition[fit_data$arm == 'CACL2_60nM'] <- '60nM CaCl2 (pH 7.5)'
fit_data$condition[fit_data$arm == 'CACL2_60mM'] <- '60mM CaCl2 (pH 7.5)'

fit_data$strn[fit_data$orf_name == 'BY4741_REF'] <- 'BY4741 (Ref.)'
fit_data$strn[fit_data$orf_name == 'BY4741'] <- 'BY4741'
fit_data$strn[fit_data$orf_name == 'VMA3'] <- '*vma3Δ::KanMX*'
fit_data$strn[fit_data$orf_name == 'VPS16'] <- '*vps16Δ::KanMX*'
fit_data$strn[fit_data$orf_name == 'YBR_CRISP'] <- '*ybr196c-aΔ*'

head(fit_data)

temp_sum <- fit_data %>%
  group_by(expt_id, density, arm, stage_id, hours, strain_id, orf_name, rep) %>%
  summarise(fit_median = median(fitness, na.rm = T), fit_mad = mad(fitness, na.rm = T),
            cs_median = median(average, na.rm = T), cs_mad = mad(average, na.rm = T),
            .groups = 'keep')
fit_data <- merge(fit_data, temp_sum, by = c('expt_id', 'density', 'arm', 'stage_id',
                                             'hours', 'strain_id', 'orf_name', 'rep'))
fit_data$fitness[fit_data$fitness < (fit_data$fit_median - 2*fit_data$fit_mad) |
                       fit_data$fitness > (fit_data$fit_median + 2*fit_data$fit_mad)] <- NA
fit_data$average[is.na(fit_data$fitness)] <- NA
fit_data <- fit_data[,c(1:14)]

fit_data.sum <- fit_data %>%
  group_by(expt_id, density, arm, condition, stage_id, hours, strain_id, orf_name, strn, rep) %>%
  summarise(fit_median = median(fitness, na.rm = T), fit_mad = mad(fitness, na.rm = T),
            cs_median = median(average, na.rm = T), cs_mad = mad(average, na.rm = T),
            .groups = 'keep')

####
fit_data$orf_name <- factor(fit_data$orf_name, levels = orf_ids)
fit_data$arm <- factor(fit_data$arm, levels = arm_ids)
fit_data$condition <- factor(fit_data$condition, levels = cnd_ids)
fit_data$strn <- factor(fit_data$strn, levels = strn_ids)

fit_data.sum$orf_name <- factor(fit_data.sum$orf_name, levels = orf_ids)
fit_data.sum$arm <- factor(fit_data.sum$arm, levels = arm_ids)
fit_data.sum$condition <- factor(fit_data.sum$condition, levels = cnd_ids)
fit_data.sum$strn <- factor(fit_data.sum$strn, levels = strn_ids)

##### PLOT
plot.sat.fit <- fit_data.sum %>%
  filter(arm != 'pH5_OLD', orf_name != 'NULL', !is.na(fit_median), hours == 65) %>%
  ggplot(aes(x = strn, y = fit_median)) +
  geom_boxplot(aes(fill = orf_name), outlier.shape = NA) +
  # stat_compare_means(method = 't.test', ref.group = 'BY4741_REF') +
  labs(x = 'Strain', y = 'Relative Fitness') +
  facet_grid(density ~ condition) +
  theme_linedraw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = titles),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 30, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.position = 'none')

ggsave("output/figs/SOLID_YBR_NACL_CACL2_SAT_FIT.jpg", plot.sat.fit,
       height = two.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)





plot.gc <- fit_data.sum %>%
  filter(arm != 'pH5_OLD', orf_name != 'NULL', !is.na(fit_median)) %>%
  ggplot() +
  stat_summary(aes(x = hours, y = cs_median, col = strn),
               fun=mean, geom="line", lwd = 1) +
  stat_summary(aes(x = hours, y = cs_median, fill = strn),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  facet_grid(density ~ condition) +
  labs(x="Hours",y="Colony Size (pix.)")+
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, margin = margin(0,0,0,0)),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = ggtext::element_markdown(size = txt),
        legend.position = c(0.88, 0.4),
        legend.spacing = unit(0.1,"mm"),
        legend.key.size = unit(5, "mm"),
        legend.background = element_rect(fill = 'transparent'))
ggsave("output/figs/SOLID_YBR_NACL_CACL2_GC.jpg", plot.gc,
       height = one.5c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)





