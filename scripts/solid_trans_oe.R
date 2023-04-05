##### YBR SOLID SCREEN RESULTS FROM TRANSLATOME EXPTS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 01/04/2023

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

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

##### GATHER DATA
data <- dbGetQuery(conn, "select * from TRANSLATOME_OE_FITNESS_DATA
                   where orf_name in ('BF_control', 'YBR196C-A')")
unique(data$condition)
data$condition <- factor(data$condition, levels = c("GA","GA2","SA","HU","HO","CF",
                                                    "DM","FL","TN"))
head(data)

ref_cs <- data %>%
  filter(data == 'cs', orf_name == 'BF_control') %>%
  group_by(attempt, data, condition, hours, plate_no, source) %>%
  summarize(ref_cs = median(average, na.rm = T), .groups = 'keep')

data <- merge(data, ref_cs, by = c('attempt', 'data', 'condition', 'hours', 'plate_no', 'source')) %>%
  mutate(norm_cs = fitness * ref_cs)

plot.gc <- data %>%
  filter(attempt != 'lin', condition != 'GA2') %>%
  ggplot() +
  stat_summary(aes(x = hours, y = norm_cs, col = orf_name),
               fun=mean, geom="line", lwd = 1) +
  stat_summary(aes(x = hours, y = norm_cs, fill = orf_name),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  facet_grid(attempt ~ condition, labeller = labeller(condition = c('GA' = 'SC-URA+GAL',
                                                                    'SA' = '1M NaCl',
                                                                    'HU' = '100 mM Hydroxyurea',
                                                                    'HO' = '3mM Hydrogen Peroxide',
                                                                    'DM' = 'DMSO',
                                                                    'FL' = '25ug/ml Fluconazole',
                                                                    'TN' = "0.6uM Tunicamycin"),
                                                      attempt = c('init' = 'Initial Screen',
                                                                  'valid' = 'Validation Screen'))) +
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
ggsave("output/figs/SOLID_YBR_TRANSLATOME_GC.jpg", plot.gc,
       height = one.5c, width = 250, units = 'mm',
       bg = 'white',
       dpi = 300)
  

