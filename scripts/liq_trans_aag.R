##### AAG-ATG MUTANT LIQUID SCREENS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 01/10/2023


##### INITIALIZE
source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/initialize.R')
platemap <- read_xlsx('/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/validation/AAG_LIQ_VAL_PLATEMAP.xlsx') %>% data.frame()
platemap$smudge[platemap$pos_id == 'A12'] <- 'Y'

platemap %>%
  group_by(orf_name, condition) %>%
  count() %>% data.frame()

##### GATHER DATA
data.aag.liq <- read_xlsx('/home/sbp29/R/Projects/adaptivefitness/rawdata/translatome/validation/AAG_LIQ_VAL_VALUES.xlsx') %>% data.frame()
data.aag.liq2 <- merge(platemap, 
                       melt(data.aag.liq, id.vars = 'Time', variable.name = 'pos_id', value.name = 'OD'),
                       by = 'pos_id', all = T)
data.aag.liq2$category[data.aag.liq2$strain_id > 100000 & data.aag.liq2$strain_id < 200000] <- 'KanMX'
data.aag.liq2$category[data.aag.liq2$strain_id > 200000 & data.aag.liq2$strain_id < 300000] <- 'AAG'
data.aag.liq2$category[data.aag.liq2$strain_id > 300000] <- 'WT'

data.aag.liq2$condition <- factor(data.aag.liq2$condition, levels = c("YPDA","SA","HU","HO","DM","FL","TN"))

data.aag.liq2$y_name <- data.aag.liq2$orf_name
data.aag.liq2$y_name[data.aag.liq2$y_name == 'orf_chr4_94133_94285_+'] <- 'YDL204W-A †'
data.aag.liq2$y_name[data.aag.liq2$y_name == 'orf_chr2_614024_614173_-'] <- 'YBR196C-A'
data.aag.liq2$y_name[data.aag.liq2$y_name == 'orf_chr14_552478_552558_-'] <- 'YNL040C-A †'

data.aag.liq2$y_name <- factor(data.aag.liq2$y_name, levels = c('YBR196C-A', 'YDL204W-A †', 'YJL077C', 'YNL040C-A †', 'YPR096C'))


##### GROWTH ANALYSIS
data.aag.gc <- SummarizeGrowthByPlate(data.aag.liq)
data.aag.gr <- NULL
for (i in 2:dim(data.aag.liq)[2]) {
  if (sum(data.aag.liq[,i] < 0.005) < 5) {
    fit0 <- fit_easylinear(data.aag.liq$Time, data.aag.liq[,i], h = 8, quota = 1);
    
    temp_res <- data.frame(sample = colnames(data.aag.liq[i]), maxgr = coef(fit0)[[3]],
                           dtime = log(2)/coef(fit0)[[3]], ltime = coef(fit0)[[4]])
    data.aag.gr <- rbind(data.aag.gr, temp_res)
  }
}

data.aag.liq.fit <- merge(platemap, merge(data.aag.gc, data.aag.gr, by = 'sample', all = T), by.x = 'pos_id', by.y = 'sample', all = T)
data.aag.liq.fit$category[data.aag.liq.fit$strain_id > 100000 & data.aag.liq.fit$strain_id < 200000] <- 'KanMX'
data.aag.liq.fit$category[data.aag.liq.fit$strain_id > 200000 & data.aag.liq.fit$strain_id < 300000] <- 'AAG'
data.aag.liq.fit$category[data.aag.liq.fit$strain_id > 300000] <- 'WT'

data.aag.liq.fit$condition <- factor(data.aag.liq.fit$condition, levels = c("YPDA","SA","HU","HO","DM","FL","TN"))

data.aag.liq.fit$y_name <- data.aag.liq.fit$orf_name
data.aag.liq.fit$y_name[data.aag.liq.fit$y_name == 'orf_chr4_94133_94285_+'] <- 'YDL204W-A †'
data.aag.liq.fit$y_name[data.aag.liq.fit$y_name == 'orf_chr2_614024_614173_-'] <- 'YBR196C-A'
data.aag.liq.fit$y_name[data.aag.liq.fit$y_name == 'orf_chr14_552478_552558_-'] <- 'YNL040C-A †'

data.aag.liq.fit$y_name <- factor(data.aag.liq.fit$y_name, levels = c('YBR196C-A', 'YDL204W-A †', 'YJL077C', 'YNL040C-A †', 'YPR096C'))


##### PLOT RESULTS
plot.od <- data.aag.liq2 %>%
  filter(is.na(smudge), strain_id > 0) %>%
  ggplot(aes(x = Time, y = OD)) +
  stat_summary(aes(col = category), geom = 'line', fun = 'median') +
  stat_summary(aes(fill = category), 
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.4) +
  scale_fill_manual(values = c('WT' = 'black',
                               'AAG' = 'purple'),
                    labels = c('WT' = 'ATG',
                               'AAG' = 'AAG')) +
  scale_color_manual(values = c('WT' = 'black',
                               'AAG' = 'purple'),
                    labels = c('WT' = 'ATG',
                               'AAG' = 'AAG')) +
  facet_grid(condition ~ y_name,
             labeller = labeller(condition = c('YPDA' = 'YPDA',
                                               'SA' = '1M NaCl',
                                               'HU' = '100 mM Hydroxyurea',
                                               'HO' = '30mM Hydrogen Peroxide',
                                               'DM' = 'DMSO',
                                               'FL' = '25ug/ml Fluconazole',
                                               'TN' = '0.5ug/ml Tunicamycin'))) +
  theme_linedraw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = txt),
        axis.text = element_text(size = txt),
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.2,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(1,0,0.1,0, "mm")))
ggsave("output/figs/LIQ_YBR_TRANS_FIT.jpg", plot.od,
       height = 250, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


