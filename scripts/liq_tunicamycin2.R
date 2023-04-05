##### YBR LIQ SCREENS WITH TUNICAMYCIN
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 01/06/2023

##### INITIALIZE
source('/home/sbp29/R/Projects/ybr/scripts/liq_initialize.R')

##### GATHER DATA
platemap <- data.frame(expt_id = '2301_06', read.csv('input/liq_screens/20230109_expt6_Platemap.csv'))
platemap$smudge[is.na(platemap$smudge)] <- 'no'
data <- data.frame(expt_id = '2301_06', read.csv('input/liq_screens/20230109_expt6_Data.csv'))

data2 <- merge(data %>% melt(id.vars = c('expt_id','Time'), variable.name = 'well_ID', value.name = 'OD'),
               platemap, by = c('expt_id','well_ID'))
unique(data2$condition)
data2$condition <- factor(data2$condition, levels = c('DMSO (pH 7.5)',"1ug/ml Tunicamycin (pH 7.5)"))

# data2$control <- 'no'
# data2$control[data2$condition != "YPD (pH 5.0)" & data2$sample %notin% c('BLANK','BY4741', '*ybr196c-aΔ*', '*vma3Δ::KanMX*')] <- 'yes'

cnd_ids <- c("DMSO (pH 7.5)", "1ug/ml Tunicamycin (pH 7.5)")
orf_ids <- c("BY4741_REF", "BY4741", "BY4742",
             "KO_PTP2_4741", "KO_PTP2_4742", "KO_HOG1_4741", "KO_HOG1_4742", "KO_VMA3_4741", "KO_VOA1_4741", "KO_VOA1_4742", "KO_VPS16_4741",
             "V1_PAM1_01", "V1_ATG_01", "V1_ATG_02", "V1_STOP_01", "V1_STOP_02", "V1_DEL_01", "V1_DEL_02",
             "V2_DEL_01", "V2_DEL_02",
             "V3_DEL_01", "V3_DEL_02", "V3_PAM1_01", "V3_PAM1_02", "V3_ATG_01", "V3_ATG_02", "V3_STOP23SYN_01", "V3_STOP23SYN_02", "V3_STOP23_01", "V3_STOP23_02",    
             "KAN_DEL_4741", "KAN_DEL_4742", "KO_DEL_4741", "KO_DEL_4742", 
             "BY4741_BOR")
data2$sample <- factor(data2$sample, levels = orf_ids)

##### GROWTH ANALYSIS
data.gc <- data.frame(expt_id = '2301_06', SummarizeGrowthByPlate(data[data$expt_id == '2301_06',-1]))
data.gr <- NULL
for (e in unique(data$expt_id)) {
  if (data$Time[data$expt_id == e][2] - data$Time[data$expt_id == e][1] > 30) {
    h <- 3
  } else {
    h <- 8
  }
  for (i in 3:dim(data[data$expt_id == e,])[2]) {
    if (sum(data[data$expt_id == e,i] < 0.005) < 5) {
      fit0 <- fit_easylinear(data$Time[data$expt_id == e], data[data$expt_id == e,i], h = h, quota = 1);
      
      temp_res <- data.frame(expt_id = e, sample = colnames(data[i]), maxgr = coef(fit0)[[3]],
                             dtime = log(2)/coef(fit0)[[3]], ltime = coef(fit0)[[4]])
      data.gr <- rbind(data.gr, temp_res)
    }
  }
}
head(data.gc)
head(data.gr)

res.liq <- merge(platemap, merge(data.gc, data.gr, by = c('expt_id','sample')), by.x = c('expt_id','well_ID'), by.y = c('expt_id','sample'))
res.liq <- res.liq %>% filter(smudge != 'yes', note == '')
res.liq <- res.liq[,colnames(res.liq) %notin% c('smudge','note')]
head(res.liq)

res.liq$condition <- factor(res.liq$condition, levels =  cnd_ids)
res.liq$sample <- factor(res.liq$sample, levels = orf_ids)

# save.image(file = 'data/ybr_liq_2301_expt6.RData')

# load(file = 'data/ybr_liq_2301_expt6.RData')
##### GROWTH CURVES
head(data2)

plot.gc <- data2 %>%
  filter(sample != 'BLANK', smudge != 'yes') %>%
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
ggsave("output/figs/LIQ_YBR_TUNI2_GC.jpg", plot.gc,
       height = one.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


##### GROWTH DYNAMICS
plot.auc <- res.liq %>%
  filter(sample != 'BLANK') %>%
  ggplot(aes(x = sample, y = auc_l)) +
  stat_summary(col = 'black',
               alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  labs(x = '', y = 'AUC') +
  facet_wrap(expt_id~condition, ncol = 1) +
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
ggsave("output/figs/LIQ_YBR_TUNI2_AUC.jpg", plot.auc,
       height = two.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


plot.dtime <- res.liq %>%
  filter(sample != 'BLANK') %>%
  ggplot(aes(x = sample, y = dtime)) +
  stat_summary(col = 'black',
               alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  scale_y_log10() +
  labs(x = '', y = 'Doubling Time (min.)') +
  facet_wrap(expt_id~condition, ncol = 1) +
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
ggsave("output/figs/LIQ_YBR_TUNI2_DTIME.jpg", plot.dtime,
       height = two.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)




##### FILTER STRAINS
data2 %>%
  filter(sample %in% c('BY4742',
                       "V1_ATG_01", "V1_ATG_02"
                       ),
         smudge != 'yes') %>%
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

res.liq %>%
  filter(sample %in% c('BY4742',
                       "V1_ATG_01", "V1_ATG_02"
  )) %>%
  ggplot(aes(x = sample, y = auc_l)) +
  stat_summary(col = 'black',
               alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  labs(x = '', y = 'AUC') +
  facet_wrap(expt_id~condition, ncol = 1) +
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


res.liq %>%
  filter(sample %in% c('BY4742',"V1_PAM1_01", "V1_ATG_01", "V1_ATG_02", "V1_STOP_01", "V1_STOP_02", "V1_DEL_01", "V1_DEL_02")) %>%
  ggplot(aes(x = sample, y = dtime)) +
  stat_summary(col = 'black',
               alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  scale_y_log10() +
  labs(x = '', y = 'Doubling Time (min.)') +
  facet_wrap(expt_id~condition, ncol = 1) +
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



