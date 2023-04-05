##### YBR LIQ SCREENS WITH NACL AND CACL2 (EXPT 38, 43 & 44)
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 12/19/2022

##### INITIALIZE
# source('/home/sbp29/R/Projects/ybr/scripts/liq_initialize.R')
# 
# ##### GATHER DATA
# platemap <- rbind(data.frame(expt_id = 38, read.csv('input/liq_screens/Expt38_Platemap.csv')),
#                   data.frame(expt_id = 43, read.csv('input/liq_screens/Expt43_Platemap.csv')),
#                   data.frame(expt_id = 44, read.csv('input/liq_screens/Expt44_Platemap.csv')))
# platemap$smudge[is.na(platemap$smudge)] <- 'no'
# data <- rbind(data.frame(expt_id = 38, read.csv('input/liq_screens/Expt38_Data.csv')),
#               data.frame(expt_id = 43, read.csv('input/liq_screens/Expt43_Data.csv')[-2]),
#              data.frame(expt_id = 44, read.csv('input/liq_screens/Expt44_Data.csv')[-2]))
# 
# data2 <- merge(data %>% melt(id.vars = c('expt_id','Time'), variable.name = 'well_ID', value.name = 'OD'),
#                platemap, by = c('expt_id','well_ID'))
# unique(data2$condition)
# data2$condition <- factor(data2$condition, levels = c('YPD (pH 5.0)',"60nM CaCl2 (pH 7.5)",'60mM CaCl2 (pH 7.5)',"1M NaCl (pH 5.0)"))
# 
# data2$control <- 'no'
# data2$control[data2$condition != "YPD (pH 5.0)" & data2$sample %notin% c('BY4741', '*ybr196c-aΔ*')] <- 'yes'
# 
# unique(data2$sample)
# data2$sample <- factor(data2$sample, levels = c("BY4741",
#                                                 "*vma3Δ::KanMX*","*vps16Δ::KanMX*","*hog1Δ::KanMX*",
#                                                 "*ybr196c-aΔ::KanMX*","*ybr196c-aΔ*",
#                                                 "BY4742","BLANK"))
# 
# ##### GROWTH ANALYSIS
# data.gc <- rbind(data.frame(expt_id = 38, SummarizeGrowthByPlate(data[data$expt_id == 38,-1])),
#                  data.frame(expt_id = 43, SummarizeGrowthByPlate(data[data$expt_id == 43,-1])),
#                  data.frame(expt_id = 44, SummarizeGrowthByPlate(data[data$expt_id == 44,-1])))
# data.gr <- NULL
# for (e in unique(data$expt_id)) {
#   if (data$Time[data$expt_id == e][2] - data$Time[data$expt_id == e][1] > 30) {
#     h <- 3
#   } else {
#     h <- 8
#   }
#   for (i in 3:dim(data[data$expt_id == e,])[2]) {
#     if (sum(data[data$expt_id == e,i] < 0.005) < 5) {
#       fit0 <- fit_easylinear(data$Time[data$expt_id == e], data[data$expt_id == e,i], h = h, quota = 1);
# 
#       temp_res <- data.frame(expt_id = e, sample = colnames(data[i]), maxgr = coef(fit0)[[3]],
#                              dtime = log(2)/coef(fit0)[[3]], ltime = coef(fit0)[[4]])
#       data.gr <- rbind(data.gr, temp_res)
#     }
#   }
# }
# head(data.gc)
# head(data.gr)
# 
# res.liq <- merge(platemap, merge(data.gc, data.gr, by = c('expt_id','sample')), by.x = c('expt_id','well_ID'), by.y = c('expt_id','sample'))
# res.liq <- res.liq %>% filter(smudge != 'yes', note == '')
# res.liq <- res.liq[,colnames(res.liq) %notin% c('smudge','note')]
# head(res.liq)
# 
# res.liq$condition <- factor(res.liq$condition, levels = c('YPD (pH 5.0)',"60nM CaCl2 (pH 7.5)",'60mM CaCl2 (pH 7.5)',"1M NaCl (pH 5.0)"))
# res.liq$sample <- factor(res.liq$sample, levels = c("BY4741",
#                                                 "*vma3Δ::KanMX*","*vps16Δ::KanMX*","*hog1Δ::KanMX*",
#                                                 "*ybr196c-aΔ::KanMX*","*ybr196c-aΔ*",
#                                                 "BY4742","BLANK"))
# 
# save.image(file = 'data/ybr_liq_expt38_43_44.RData')

load(file = 'data/ybr_liq_expt38_43_44.RData')
##### GROWTH CURVES
head(data2)

plot.gc <- data2 %>%
  filter(sample != 'BLANK', smudge != 'yes', sample != 'BY4742') %>%
  ggplot() +
  stat_summary(aes(x = Time, y = OD, col = sample),
               fun=mean, geom="line", lwd = 1) +
  stat_summary(aes(x = Time, y = OD, fill = sample),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  facet_grid(expt_id~condition) +
  labs(x="Minute",y="OD")+
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, margin = margin(0,0,0,0)),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = ggtext::element_markdown(size = txt),
        legend.position = 'bottom',
        legend.spacing = unit(0.1,"mm"),
        legend.key.size = unit(5, "mm"),
        legend.background = element_rect(fill = 'transparent'))
ggsave("output/figs/LIQ_YBR_NACL_CACL2_GC.jpg", plot.gc,
       height = two.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


##### GROWTH DYNAMICS
plot.auc <- res.liq %>%
  filter(sample != 'BLANK', sample != 'BY4742') %>%
  ggplot(aes(x = sample, y = auc_l)) +
  stat_summary(col = 'black',
               alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  labs(x = '', y = 'AUC') +
  facet_grid(expt_id~condition) +
  theme_linedraw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 30, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.2,"mm"))
ggsave("output/figs/LIQ_YBR_NACL_CACL2_AUC.jpg", plot.auc,
       height = two.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


plot.dtime <- res.liq %>%
  filter(sample != 'BLANK', sample != 'BY4742') %>%
  ggplot(aes(x = sample, y = dtime)) +
  stat_summary(col = 'black',
               alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  scale_y_log10() +
  labs(x = '', y = 'Doubling Time (min.)') +
  facet_grid(expt_id~condition) +
  theme_linedraw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 30, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.2,"mm"))
ggsave("output/figs/LIQ_YBR_NACL_CACL2_DTIME.jpg", plot.dtime,
       height = two.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


