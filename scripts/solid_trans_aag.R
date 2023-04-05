##### AAG-ATG MUTANT ANALYSIS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 01/09/2023

##### INITIALIZE
source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/initialize.R')
source('/home/sbp29/R/Projects/adaptivefitness/scripts/translatome/fit_sum.R')

##### GATHER DATA
data.aag <- dbGetQuery(conn, 'select * from TR_AAG_CS_FIT')
head(data.aag)
unique(data.aag$attempt)
data.aag$attempt <- factor(data.aag$attempt, levels = c('ONE','TWO'))
unique(data.aag$condition)
data.aag$condition <- factor(data.aag$condition, levels = c("YPDA","SA","HU","HO",
                                                            "DM","FL","TN"))
data.aag <- data.aag %>%
  mutate(id = paste(attempt, condition, orf_name, pos, rep, sep = '_'))
data.aag$hours <- as.numeric(data.aag$hours)

data.aag$y_name <- data.aag$orf_name
unique(data.aag$orf_name)
data.aag$y_name[data.aag$y_name == 'orf_chr4_94133_94285_+'] <- 'YDL204W-A †'
data.aag$y_name[data.aag$y_name == 'orf_chr2_614024_614173_-'] <- 'YBR196C-A'
data.aag$y_name[data.aag$y_name == 'orf_chr4_593209_593304_-'] <- 'YDR073C-A †'
data.aag$y_name[data.aag$y_name == 'orf_chr7_523246_523353_-'] <- 'YGR016C-A †'
data.aag$y_name[data.aag$y_name == 'orf_chr14_552478_552558_-'] <- 'YNL040C-A †'
unique(data.aag$y_name)
data.aag$y_name <- factor(data.aag$y_name, levels = c('YBR196C-A', 'YDL204W-A †', 'YDR073C-A †', 'YGR016C-A †', 'YJL077C', 'YNL040C-A †', 'YPR096C',
                                                      'FY4', 'BY4741', 'BY4741_HO'))

unique(data.aag$category)
data.aag$category <- factor(data.aag$category, levels = c('WT','AAG','KanMX'))


##### GROWTH CURVES
data.aag %>%
  filter(n > 4, orf_name %in% c('orf_chr2_614024_614173_-'), category %in% c('WT','AAG')) %>%
  ggplot(aes(x = hours, y = average, col = category, fill = category)) +
  stat_summary(fun=mean, geom="line", lwd = 1) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  facet_wrap(attempt ~ condition, nrow = 2)


data.aag %>%
  filter(n > 4, orf_name %in% c('orf_chr2_614024_614173_-'), category %in% c('WT','AAG'), hours == 90) %>%
  group_by(attempt, condition, category) %>% count() %>% data.frame()

stats.fit <- compare_means(fitness ~ category, data.aag %>%
                               filter(n > 4, category %in% c('WT','AAG'), hours == 90),
                             group.by = c('condition','orf_name','y_name'), method = 'kruskal')

plot.fit <- data.aag %>%
  filter(n > 4,
         category %in% c('WT','AAG'),
         hours == 90) %>%
  ggplot(aes(x = y_name, y = fitness)) +
  # geom_jitter(size = 0.2) +
  geom_violin(aes(fill = category), alpha = 0.5, draw_quantiles = c(0.25, 0.5, 0.75),
              size = 0.2) + 
  geom_text(data = stats.fit, aes(x = y_name, y = 1.5, label = p.signif), size = 2) +
  # stat_compare_means(method = 'kruskal',
  #                    aes(label = paste0("p = ", ..p.format..)), size = 2, hjust = 0) +
  scale_fill_manual(values = c('WT' = 'black',
                               'AAG' = 'purple'),
                    labels = c('WT' = 'ATG',
                               'AAG' = 'AAG')) +
  labs(x = '', y = 'Fitness') +
  facet_wrap(.~condition, ncol = 1,
             labeller = labeller(condition = c('YPDA' = 'YPDA',
                                               'SA' = '1M NaCl',
                                               'HU' = '100 mM Hydroxyurea',
                                               'HO' = '30mM Hydrogen Peroxide',
                                               'DM' = 'DMSO',
                                               'FL' = '25ug/ml Fluconazole',
                                               'TN' = "0.5ug/ml Tunicamycin"))) +
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
ggsave("output/figs/SOLID_YBR_TRANS_FIT.jpg", plot.fit,
       height = 220, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)



##### AAG - ATG DIFFERENCE
data.aag.diff <- NULL
for (a in unique(data.aag$attempt)) {
  for (c in unique(data.aag$condition[data.aag$attempt == a])) {
    for (y in unique(data.aag$y_name[data.aag$attempt == a & data.aag$condition == c])) {
      temp.fit.aag <- data.aag$fitness[data.aag$attempt == a & data.aag$condition == c & 
                                         data.aag$y_name == y & data.aag$hours == 90 &
                                         data.aag$category == 'AAG']
      temp.fit.atg <- data.aag$fitness[data.aag$attempt == a & data.aag$condition == c & 
                                         data.aag$y_name == y & data.aag$hours == 90 &
                                         data.aag$category == 'WT']
      
      n <- length(temp.fit.aag[!is.na(temp.fit.aag)]) + length(temp.fit.atg[!is.na(temp.fit.atg)])/2
      
      temp.fit.diff <- c(sapply(temp.fit.aag,"-",temp.fit.atg))
      
      data.aag.diff <- rbind(data.aag.diff,
                             data.frame(attempt = a, condition = c,
                                        orf_name = unique(data.aag$orf_name[data.aag$attempt == a & data.aag$condition == c &data.aag$y_name == y]),
                                        y_name = y,
                                        mean_fit_diff = mean(temp.fit.diff, na.rm = T),
                                        median_fit_diff = median(temp.fit.diff, na.rm = T),
                                        std_fit_diff = sd(temp.fit.diff, na.rm = T),
                                        se_fit_diff = sd(temp.fit.diff, na.rm = T)/sqrt(n)))
    }
  }
}
head(data.aag.diff)
data.aag.diff <- data.aag.diff %>%
  group_by(condition, orf_name, y_name) %>%
  summarise(mean_fit_diff = mean(mean_fit_diff),
            median_fit_diff = mean(median_fit_diff),
            std_fit_diff = mean(std_fit_diff),
            se_fit_diff = mean(se_fit_diff),
            .groups = 'keep') %>%
  data.frame()

stats.fit <- merge(stats.fit, data.aag.diff, by = c('condition','orf_name','y_name'))


plot.fit.diff <- stats.fit %>%
  ggplot(aes(x = condition, y = mean_fit_diff)) +
  stat_summary(col = 'black',
               alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  geom_errorbar(aes(x=condition,
                    ymin=mean_fit_diff-(1.96*se_fit_diff),
                    ymax=mean_fit_diff+(1.96*se_fit_diff)),
                width=0.4, colour="black", alpha=0.9, size=0.5) +
  geom_text(aes(x = condition, y = 0.3, label = p.signif), size = 2) +
  scale_x_discrete(limits = c('YPDA', 'SA', 'HU', 'HO',
                              'DM','FL','TN'),
                   labels = c('YPDA' = 'YPDA',
                              'SA' = '1M NaCl',
                              'HU' = '100 mM Hydroxyurea',
                              'HO' = '30mM Hydrogen Peroxide',
                              'DM' = 'DMSO',
                              'FL' = '25ug/ml Fluconazole',
                              'TN' = "0.5ug/ml Tunicamycin")) +
  labs(x = '', y = 'Fitness difference between AAG and ATG mutant') +
  facet_wrap(.~y_name, ncol = 2) +
  theme_linedraw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = txt),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.margin = margin(0,0,0,0),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.2,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(1,0,0.1,0, "mm")))

ggsave("output/figs/SOLID_YBR_TRANS_FIT_DIFF.jpg", plot.fit.diff,
       height = 220, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)

