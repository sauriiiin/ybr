##### YBR196C-A Protein Western Blot Quantifications
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 02/09/2023

##### INITIALIZE
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(dplyr)
library(stringr)
library(reshape2)
`%notin%` <- Negate(`%in%`)

##### LOAD DATA
doa10_steady_state <- read.csv(file = 'input/ally_wb_quants/doa10_steady_state.csv')
chx_treatment <- read.csv(file = 'input/ally_wb_quants/protein_synthesis_inhibition.csv')
colnames(chx_treatment) <- c("Background", "Minutes", "ProteinAmount", "Sample", "Media")
proteosome_inhibition <- read.csv(file = 'input/ally_wb_quants/proteosome_inhibition.csv')

background_id <- c("WT","*hrd1Δ*","*doa10Δ*","*doa10Δ hrd1Δ*","*pdr5Δ*")
sample_id <- c("YBR196C-A-mNG", "mNG", "WT", "*doa10Δ*")

doa10_steady_state$Sample <- factor(doa10_steady_state$Sample, level = sample_id)
chx_treatment$Sample <- factor(chx_treatment$Sample, level = sample_id)
proteosome_inhibition$Sample <- factor(proteosome_inhibition$Sample, level = sample_id)

doa10_steady_state$Background <- str_replace_all(doa10_steady_state$Background, 'D', 'Δ')
doa10_steady_state$Background <- factor(doa10_steady_state$Background, level = background_id)
chx_treatment$Background <- str_replace_all(chx_treatment$Background, 'D', 'Δ')
chx_treatment$Background <- factor(chx_treatment$Background, level = background_id)
proteosome_inhibition$Background <- str_replace_all(proteosome_inhibition$Background, 'D', 'Δ')
proteosome_inhibition$Background <- factor(proteosome_inhibition$Background, level = background_id)



doa10_steady_state <- merge(doa10_steady_state, doa10_steady_state %>%
  filter(Background == 'WT') %>%
  group_by(Sample) %>%
  summarise(.groups = 'keep',
            WT_ProteinAmount = median(ProteinAmount)),
  by = 'Sample') %>%
  mutate(RelativeProteinAbundance = ProteinAmount/WT_ProteinAmount)

plot.doa10.steady.state <- doa10_steady_state %>%
  ggplot(aes(x = Background, y = RelativeProteinAbundance)) +
  stat_summary(col = 'black', fill = '#BDBDBD',
               alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  stat_compare_means(method = 't.test', ref.group = 'WT',
                     aes(label = paste0("p = ", ..p.format..))) +
  labs(y = 'Relative Protein Abundance') +
  facet_wrap(.~Sample) +
  theme_linedraw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = ggtext::element_markdown(size = 8),
        axis.text.y = element_text(size = 8))
ggsave("output/figs/ybr_presentation/protein_abundance_doa10_steady_state.jpg", plot.doa10.steady.state,
       height = 100, width = 200, units = 'mm',
       bg = 'white',
       dpi = 300)


plot.chx.treatment <- chx_treatment %>%
  ggplot(aes(x = Minutes, y = ProteinAmount)) +
  # stat_summary(aes(group = Background, fill = Background),
  #              fun.data=mean_se, fun.args = list(mult=1), geom="ribbon", alpha = 0.4) +
  stat_summary(aes(group = Background, col = Background), fun=mean, geom="line", lwd = 0.7) +
  stat_summary(aes(group = Background, col = Background), fun=mean, geom="point", size = 2) +
  stat_summary(aes(group = Background, col = Background),
               fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1),
               size = 1) +
  labs(x = 'Time after CHX (min.)', y = 'Relative Protein Abundance') +
  theme_linedraw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.position = 'bottom',
        legend.title = element_text(size = 10),
        legend.text = ggtext::element_markdown(size = 8),
        legend.margin = margin(-10,0,0,0))
ggsave("output/figs/ybr_presentation/protein_abundance_chx_treatment.jpg", plot.chx.treatment,
       height = 100, width = 200, units = 'mm',
       bg = 'white',
       dpi = 300)


plot.proteosome.inhibition <- proteosome_inhibition %>%
  ggplot(aes(x = Minutes, y = ProteinAmount)) +
  stat_summary(aes(group = Media, col = Media), fun=mean, geom="line", lwd = 0.7) +
  stat_summary(aes(group = Media, col = Media), fun=mean, geom="point", size = 2) +
  stat_summary(aes(group = Media, col = Media),
               fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1),
               size = 1) +
  labs(x = 'Time after CHX (min.)', y = 'Relative Protein Abundance in *pdr5Δ*') +
  theme_linedraw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.title.y = ggtext::element_markdown(size = 10),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.margin = margin(-10,0,0,0))
ggsave("output/figs/ybr_presentation/protein_abundance_proteosome_inhibition.jpg", plot.proteosome.inhibition,
       height = 100, width = 200, units = 'mm',
       bg = 'white',
       dpi = 300)
