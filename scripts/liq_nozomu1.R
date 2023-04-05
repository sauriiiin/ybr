##### YBR LIQ SCREENS WITH CaCl2, NaCl and Tunicamycin
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 01/09/2023

##### INITIALIZE
source('/home/sbp29/R/Projects/ybr/scripts/liq_initialize.R')

##### GATHER DATA
platemap <- rbind(data.frame(expt_id = '2301_11', read.csv('input/liq_screens/20230124_expt11_Platemap.csv')))
platemap$smudge[is.na(platemap$smudge)] <- 'no'

data <- rbind(data.frame(expt_id = '2301_11', read.csv('input/liq_screens/20230124_expt11_Data.csv')))

data2 <- merge(data %>% melt(id.vars = c('expt_id','Time'), variable.name = 'well_ID', value.name = 'OD'),
               platemap, by = c('expt_id','well_ID'))
unique(data2$condition)
cnd_ids <- c("YPDA (pH 6.8)",
             "1M NaCl (pH 6.8)",
             "1M NaCl + 5mM CaCl2 (pH 6.8)",
             "1M NaCl + 60mM CaCl2 (pH 6.8)")
data2$condition <- factor(data2$condition, levels = cnd_ids)

data2$control <- NA

unique(data2$sample)
orf_ids <-  c("empty","YBROE",
              "YBRDel_empty","YBRDel_YBROE",
              "BLANK")
data2$sample <- factor(data2$sample, levels = orf_ids)

##### GROWTH ANALYSIS
data.gc <- rbind(data.frame(expt_id = '2301_11', SummarizeGrowthByPlate(data[data$expt_id == '2301_11',-1])))
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

res.liq$condition <- factor(res.liq$condition, levels = cnd_ids)
res.liq$sample <- factor(res.liq$sample, levels = orf_ids)

# save.image(file = 'data/ybr_liq_expt44_03_05.RData')
# load(file = 'data/ybr_liq_expt44_03_05.RData')
##### GROWTH CURVES
head(data2)

plot.gc <- data2 %>%
  filter(sample != 'BLANK', smudge != 'yes') %>%
  ggplot() +
  stat_summary(aes(x = Time, y = OD, col = sample),
               fun=mean, geom="line", lwd = 1) +
  stat_summary(aes(x = Time, y = OD, fill = sample),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  facet_wrap(expt_id~condition,dir = 'v',nrow = 2) +
  labs(x="Minute",y="OD")+
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, margin = margin(0,0,0,0)),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        legend.text = ggtext::element_markdown(size = txt),
        legend.position = 'bottom',
        legend.spacing = unit(0.1,"mm"),
        legend.key.size = unit(5, "mm"),
        legend.background = element_rect(fill = 'transparent'))
ggsave("output/figs/LIQ_YBR_NOZOMU_GC.jpg", plot.gc,
       height = two.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


##### GROWTH DYNAMICS
# stats.auc_l <- compare_means(auc_l ~ condition, res.liq %>%
#                                filter(sample != 'BLANK'),
#                              group.by = c('expt_id','sample'), ref.group = "YPDA (pH 6.8)", method = 't.test')
# 
# stats.dtime <- compare_means(dtime ~ condition, res.liq %>%
#                                filter(sample != 'BLANK'),
#                              group.by = c('expt_id','sample'), ref.group = "YPDA (pH 6.8)", method = 't.test')

sum.liq <- res.liq %>%
  filter(sample != 'BLANK') %>%
  group_by(expt_id, condition, sample) %>%
  summarise(dtime = median(dtime, na.rm = T),
            auc_l = median(auc_l, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

res.liq <- merge(res.liq, sum.liq %>% filter(sample == 'empty'),
      by = c('expt_id', 'condition'), suffixes = c('','_control')) %>%
  mutate(norm_dtime = dtime/dtime_control, norm_auc = auc_l/auc_l_control)

plot.auc <- res.liq %>%
  filter(sample != 'BLANK') %>%
  ggplot(aes(x = condition, y = norm_auc)) +
  stat_summary(col = 'black',
               alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  # geom_text(data = stats.auc_l, aes(x = group2, y = 1.2, label = p.signif)) +
  labs(x = '', y = 'AUC') +
  facet_wrap(expt_id~sample,dir = 'v',nrow = 1, scales = 'free_x') +
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
ggsave("output/figs/LIQ_YBR_NOZOMU_AUC.jpg", plot.auc,
       height = one.5c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


plot.dtime <- res.liq %>%
  filter(sample != 'BLANK') %>%
  ggplot(aes(x = condition, y = dtime)) +
  stat_summary(col = 'black',
               alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  # geom_text(data = stats.dtime, aes(x = group2, y = 1300, label = p.signif)) +
  # scale_y_log10() +
  labs(x = '', y = 'Doubling Time (min.)') +
  facet_wrap(expt_id~sample,dir = 'v',nrow = 1, scales = 'free_x') +
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
ggsave("output/figs/LIQ_YBR_NOZOMU_DTIME.jpg", plot.dtime,
       height = two.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)



plot.res <- plot_grid(plot.gc + facet_wrap(expt_id~condition, nrow = 1),
                      plot.auc + facet_wrap(expt_id~sample, nrow = 1),
                      plot.dtime + facet_wrap(expt_id~sample, nrow = 1),
                      nrow = 3,
                      labels = c('A','B','C'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave("output/figs/LIQ_YBR_NOZOMU_RES.jpg", plot.res,
       height = 240, width = 280, units = 'mm',
       bg = 'white',
       dpi = 300)

