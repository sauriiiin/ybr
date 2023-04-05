##### YBR SOLID SCREENS WITH NACL, CACL2 and TUNI (EXPT 1)
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 01/17/2022

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

file_info <- readxl::read_xlsx(path = '/home/sbp29/RAW_Data/YBR_NACL_CACL_TUNI/20230106_YBRs_INFO.xlsx') %>% data.frame()
head(file_info)

arm_ids <- c("YPDA_5", "NACL_5", "YPDA_75", "CACL2_75", "DMSO_75", "TUNI_75")
cnd_ids <- c("YPDA pH5.0", "1M NaCl pH 5.0", "YPDA pH7.5", "60mM CaCl2 pH 7.5", "DMSO pH 7.5", "1ug/ml Tunicamycin pH 7.5")
orf_ids <- c("BY4741_REF", "BY4741", "BY4742",
             "KO_PTP2_4741", "KO_PTP2_4742", "KO_HOG1_4741", "KO_HOG1_4742", "KO_VMA3_4741", "KO_VOA1_4741", "KO_VOA1_4742", "KO_VPS16_4741",
             "V1_PAM1_01", "V1_ATG_01", "V1_ATG_02", "V1_STOP_01", "V1_STOP_02", "V1_DEL_01", "V1_DEL_02",
             "V2_DEL_01", "V2_DEL_02",
             "V3_DEL_01", "V3_DEL_02", "V3_PAM1_01", "V3_PAM1_02", "V3_ATG_01", "V3_ATG_02", "V3_STOP23SYN_01", "V3_STOP23SYN_02", "V3_STOP23_01", "V3_STOP23_02",    
             "KAN_DEL_4741", "KAN_DEL_4742", "KO_DEL_4741", "KO_DEL_4742", 
             "BY4741_BOR")

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

fit_data$condition[fit_data$arm == "YPDA_5"] <- "YPDA pH5.0"
fit_data$condition[fit_data$arm == "NACL_5"] <- "1M NaCl pH 5.0"
fit_data$condition[fit_data$arm == "YPDA_75"] <- "YPDA pH7.5"
fit_data$condition[fit_data$arm == "CACL2_75"] <- "60mM CaCl2 pH 7.5"
fit_data$condition[fit_data$arm == "DMSO_75"] <- "DMSO pH 7.5"
fit_data$condition[fit_data$arm == "TUNI_75"] <- "1ug/ml Tunicamycin pH 7.5"

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
fit_data <- fit_data[,c(1:13)]

fit_data.sum <- fit_data %>%
  group_by(expt_id, density, arm, condition, stage_id, hours, strain_id, orf_name, rep) %>%
  summarise(fit_median = median(fitness, na.rm = T), fit_mad = mad(fitness, na.rm = T),
            cs_median = median(average, na.rm = T), cs_mad = mad(average, na.rm = T),
            .groups = 'keep')

####
fit_data$orf_name <- factor(fit_data$orf_name, levels = orf_ids)
fit_data$arm <- factor(fit_data$arm, levels = arm_ids)
fit_data$condition <- factor(fit_data$condition, levels = cnd_ids)

fit_data.sum$orf_name <- factor(fit_data.sum$orf_name, levels = orf_ids)
fit_data.sum$arm <- factor(fit_data.sum$arm, levels = arm_ids)
fit_data.sum$condition <- factor(fit_data.sum$condition, levels = cnd_ids)

##### PLOT
plot.sat.fit <- fit_data %>%
  filter(!is.na(fitness), hours == 76) %>%
  ggplot(aes(x = orf_name, y = fitness)) +
  geom_boxplot(aes(fill = orf_name), outlier.shape = NA) +
  labs(x = 'Strain', y = 'Relative Fitness') +
  facet_wrap(.~ condition, ncol = 1) +
  theme_linedraw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = titles),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 30, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.position = 'none')

ggsave("output/figs/SOLID_YBR_NACL_CACL2_TUNI_SAT_FIT.jpg", plot.sat.fit,
       height = two.c*1.5, width = two.c*1.5, units = 'mm',
       bg = 'white',
       dpi = 300)





plot.gc <- fit_data.sum %>%
  filter(!is.na(fit_median)) %>%
  ggplot() +
  stat_summary(aes(x = hours, y = cs_median, col = orf_name),
               fun=mean, geom="line", lwd = 1) +
  stat_summary(aes(x = hours, y = cs_median, fill = orf_name),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  facet_wrap(.~condition, ncol =1) +
  labs(x="Hours",y="Colony Size (pix.)")+
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, margin = margin(0,0,0,0)),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = ggtext::element_markdown(size = txt),
        legend.position = 'right',
        legend.spacing = unit(0.1,"mm"),
        legend.key.size = unit(5, "mm"),
        legend.background = element_rect(fill = 'transparent'))
ggsave("output/figs/SOLID_YBR_NACL_CACL2_TUNI_GC.jpg", plot.gc,
       height = two.c*1.5, width = two.c*1.5, units = 'mm',
       bg = 'white',
       dpi = 300)





