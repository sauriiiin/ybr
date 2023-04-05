##### YBR LIQ SCREENS WITH NACL AND CACL2
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 12/07/2022

##### INITIALIZE
source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/initialize.R')

##### GATHER DATA
platemap <- read.csv('input/liq_screens/ybr_nacl_cacl2_platemap.csv')

data <- read.csv('input/liq_screens/YBR_NaCl_CaCl2.csv')
data2 <- merge(data %>% melt(id.vars = c('Time'), variable.name = 'well_ID', value.name = 'OD'),
               platemap, by = 'well_ID')
unique(data2$condition)
data2$condition[data2$condition == '60nM CaCl2 (pH 7.5)'] <- '60mM CaCl2 (pH 7.5)'
data2$condition <- factor(data2$condition, levels = c('YPD (pH 5.0)', '1M NaCl (pH 5.0)', '60mM CaCl2 (pH 7.5)'))

data2$control[data2$condition == 'YPD (pH 5.0)'] <- ''
data2$control[data2$control == ''] <- 'no'

data2$background[str_detect(data2$sample,'BY4741')] <- 'BY4741'
data2$background[str_detect(data2$sample,'BY4742')] <- 'BY4742'

data2$sample <- str_remove(data2$sample, 'BY4741-')
data2$sample <- str_remove(data2$sample, 'BY4742-')
unique(data2$sample)

data2$background <- factor(data2$background, levels = c('BY4741','BY4742'))
data2$sample <- factor(data2$sample, levels = c("BY4741","BY4742",
                                                "*hog1Δ::KanMX*","*vps16Δ::KanMX*","*vma3Δ::KanMX*",
                                                "*ybr196c-aΔ::KanMX*","*ybr196c-aΔ*",
                                                "BLANK"))

##### GROWTH ANALYSIS
# data.gc <- SummarizeGrowthByPlate(data)
# data.gr <- NULL
# for (i in 2:dim(data)[2]) {
#   if (sum(data[,i] < 0.005) < 5) {
#     fit0 <- fit_easylinear(data$Time, data[,i], h = 8, quota = 1);
#     
#     temp_res <- data.frame(sample = colnames(data[i]), maxgr = coef(fit0)[[3]],
#                            dtime = log(2)/coef(fit0)[[3]], ltime = coef(fit0)[[4]])
#     data.gr <- rbind(data.gr, temp_res)
#   }
# }
# head(data.gc)
# head(data.gr)
# res.liq <- merge(platemap, merge(data.gc, data.gr, by = 'sample'), by.x = 'well_ID', by.y = 'sample')
# res.liq <- res.liq %>% filter(smudge != 'yes', note == '')
# res.liq <- res.liq[,colnames(res.liq) %notin% c('smudge','note')]

##### GROWTH CURVES
head(data2)

data2 %>%
  filter(sample != 'BLANK', smudge != 'yes') %>%
  ggplot() +
  stat_summary(aes(x = Time, y = OD, col = sample),
               fun=mean, geom="line", lwd = 1) +
  stat_summary(aes(x = Time, y = OD, fill = sample),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  facet_wrap(.~condition*control) +
  labs(x="Minute",y="OD")+
  # scale_color_manual(values=c("purple","black")) +
  # scale_fill_manual(values=c("purple","black")) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, margin = margin(0,0,0,0)),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = ggtext::element_markdown(size = txt),
        legend.position = c(0.9, 0.25),
        legend.spacing = unit(0.1,"mm"),
        legend.key.size = unit(5, "mm"),
        legend.background = element_rect(fill = 'transparent'))


data2 %>%
  filter((sample != 'BLANK' & smudge != 'yes') &
           (control == 'yes' | sample == 'BY4741')) %>%
  ggplot() +
  stat_summary(aes(x = Time, y = OD, col = sample),
               fun=mean, geom="line", lwd = 1) +
  stat_summary(aes(x = Time, y = OD, fill = sample),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  facet_wrap(.~condition) +
  labs(x="Minute",y="OD")+
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, margin = margin(0,0,0,0)),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = ggtext::element_markdown(size = txt),
        legend.position = c(0.9, 0.25),
        legend.spacing = unit(0.1,"mm"),
        legend.key.size = unit(5, "mm"),
        legend.background = element_rect(fill = 'transparent'))


data2 %>%
  filter((Time == max(data2$Time) & sample != 'BLANK' & smudge != 'yes') &
                       (control == 'yes' | sample == 'BY4741')) %>%
  ggplot(aes(x = sample, y = OD)) +
  geom_boxplot() +
  facet_wrap(.~condition) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 30, hjust = 1),
        axis.text.y = element_text(size = txt)) +
  coord_cartesian(ylim = c(0.8,2))


data2 %>%
  filter((Time == max(data2$Time) & sample != 'BLANK' & smudge != 'yes') &
           ((str_detect(sample, 'ybr') & !str_detect(sample, 'Kan')) | sample %in% c('BY4741','BY4742'))) %>%
  ggplot() +
  stat_summary(aes(x = Time, y = OD, col = sample),
               fun=mean, geom="line", lwd = 1) +
  stat_summary(aes(x = Time, y = OD, fill = sample),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  facet_wrap(.~condition) +
  labs(x="Minute",y="OD")+
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, margin = margin(0,0,0,0)),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = ggtext::element_markdown(size = txt),
        legend.position = c(0.9, 0.25),
        legend.spacing = unit(0.1,"mm"),
        legend.key.size = unit(5, "mm"),
        legend.background = element_rect(fill = 'transparent'))



data2 %>%
  filter((Time == max(data2$Time) & sample != 'BLANK' & smudge != 'yes') &
           ((str_detect(sample, 'ybr') & !str_detect(sample, 'Kan')) | sample %in% c('BY4741','BY4742'))) %>%
  ggplot(aes(x = sample, y = OD)) +
  stat_summary(col = 'black', fill = '#D6E4F4',
               alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  # stat_compare_means(method = 't.test', label = 'p.format',
  #                    # label.x = 1.5, label.y = 6.5,
  #                    hjust = 0.5, size = 2.5,
  #                    ref.group = 'BY4741') +
  labs(x = '', y = 'OD') +
  facet_wrap(.~condition) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 30, hjust = 1),
        axis.text.y = element_text(size = txt)) +
  coord_cartesian(ylim = c(0.8,2))



plot.sat.od <- data2 %>%
  # filter(!(str_detect(sample, 'ybr') & str_detect(sample, 'Kan'))) %>%
  filter(Time == max(data2$Time) & sample != 'BLANK' & smudge != 'yes' & background == 'BY4741') %>%
  ggplot(aes(x = sample, y = OD)) +
  stat_summary(aes(fill = control), col = 'black',
               alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  # stat_compare_means(method = 't.test', label = 'p.format',
  #                    # label.x = 1.5, label.y = 6.5,
  #                    hjust = 0.5, size = 2.5,
  #                    ref.group = 'BY4741') +
  labs(x = '', y = 'OD<sub>600</sub> at Saturation') +
  scale_fill_manual(values = c('no' = '#3F51B5',
                               'yes' = '#FFC107'),
                    labels = c('no' = 'Sample',
                               'yes' = 'Control')) +
  facet_wrap(.~condition) +
  theme_linedraw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 30, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.2,"mm")) +
  coord_cartesian(ylim = c(0.8,2))

ggsave("output/figs/LIQ_YBR_NACL_CACL2_SAT_OD.jpg", plot.sat.od,
       height = one.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)
