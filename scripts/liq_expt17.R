##### YBR LIQ SCREENS (EXPT 17)
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 01/04/2023

##### INITIALIZE
source('/home/sbp29/R/Projects/ybr/scripts/liq_initialize.R')

##### GATHER DATA
platemap <- data.frame(expt_id = 17, read.csv('input/liq_screens/Expt17_Platemap.csv'))
platemap$sample[platemap$sample == "*ybr196c-a::STOP*"] <- '*ybrSTOP*'
data <- data.frame(expt_id = 17, read.csv('input/liq_screens/Expt17_Data.csv'))

data2 <- merge(data %>% melt(id.vars = c('expt_id','Time'), variable.name = 'well_ID', value.name = 'OD'),
               platemap, by = c('expt_id','well_ID'))
unique(data2$condition)
data2$condition <- factor(data2$condition, levels = c('YPD',"1M NaCl","30 mM H2O2",'DMSO',"25ug/ml Fluconazole","0.6uM Tunicamycin"))

unique(data2$sample)
data2$sample <- factor(data2$sample, levels = c("BY4741",
                                                "*ybr196c-aΔ::KanMX*","*ybr196c-aΔ*","*ybrSTOP*",
                                                "BLANK"))

##### GROWTH ANALYSIS
data.gc <- data.frame(expt_id = 17, SummarizeGrowthByPlate(data[data$expt_id == 17,-1]))
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

res.liq$condition <- factor(res.liq$condition, levels = c('YPD',"1M NaCl","30 mM H2O2",'DMSO',"25ug/ml Fluconazole","0.6uM Tunicamycin"))
res.liq$sample <- factor(res.liq$sample, levels = c("BY4741",
                                                "*ybr196c-aΔ::KanMX*","*ybr196c-aΔ*","*ybrSTOP*",
                                                "BLANK"))

save.image(file = 'data/ybr_liq_expt17.RData')


load(file = 'data/ybr_liq_expt17.RData')
##### GROWTH CURVES
head(data2)

# plot.gc <-
data2 %>%
  filter(sample != 'BLANK', smudge != 'yes') %>%
  ggplot() +
  stat_summary(aes(x = Time, y = OD, col = sample, linetype = rep_ID),
               fun=mean, geom="line", lwd = 1) +
  stat_summary(aes(x = Time, y = OD, fill = sample),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  facet_wrap(.~condition, nrow = 2) +
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
ggsave("output/figs/LIQ_YBR_EXPT17_GC.jpg", plot.gc,
       height = two.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


##### GROWTH DYNAMICS
head(res.liq)
plot.auc <- res.liq %>%
  filter(sample != 'BLANK') %>%
  ggplot(aes(x = sample, y = auc_l)) +
  stat_summary(col = 'black',
               alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  labs(x = '', y = 'AUC') +
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
        legend.box.spacing = unit(0.2,"mm")) #+
  # coord_cartesian(ylim = c(0.8,2))
ggsave("output/figs/LIQ_YBR_EXPT17_AUC.jpg", plot.auc,
       height = two.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)

