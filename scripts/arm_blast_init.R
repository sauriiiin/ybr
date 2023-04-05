
library(dplyr)
library(ggplot2)

arm_domains <- read.csv(file = 'input/ms_data/HA_36/arm_proteins_all_domains.csv',
                        col.names = c('orf_name','std_name','org_name','domain_name','domain_description','start','end'))
arm_domains <- arm_domains %>% filter(domain_name == 'SSF48371')

arm_blast <- read.csv(file = 'input/ms_data/HA_36/ARM_proteins_blast.csv')
arm_blast <- arm_blast %>% mutate(id = seq(1,103))
arm_blast$id2 <- arm_blast$id
head(arm_blast)

for (q in unique(arm_blast$query)) {
  for (s in unique(arm_blast$subject[arm_blast$query == q])) {
    arm_blast$id2[arm_blast$query == q & arm_blast$subject == s] <- seq(1,length(arm_blast$id[arm_blast$query == q & arm_blast$subject == s]))
  }
}

for (q in unique(arm_blast$query)) {
  q_plot <- arm_blast %>%
    filter(query == q) %>%
    ggplot() +
    geom_segment(aes(x = 0, xend = query_length, y = '', yend = ''), lwd = 1) +
    geom_segment(aes(x = query_start, xend = query_end,
                     y = as.character(id), yend = as.character(id),
                     col = subject), lwd = 1) +
    geom_text(aes(x = (query_start+query_end)/2, y = as.character(id), label = sprintf('%0.3f', e_value)),
              size = 2.5, nudge_y = 0.6) +
    geom_segment(data = arm_domains %>% filter(std_name == q) %>% mutate(query = std_name),
                 aes(x = start, xend = end,
                     y = '999', yend = '999'),col = 'grey', lwd = 1) +
    scale_x_continuous(breaks = seq(0,5000,500)) +
    labs(x = 'Query AA Seq') +
    facet_wrap(.~query) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 11),
          axis.text.x = element_text(size = 8),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position = 'None')
  
  s_plot <- arm_blast %>%
    filter(query == q) %>%
    ggplot() +
    geom_segment(aes(x = 0, xend = subject_length, y = 0, yend = 0), lwd = 1) +
    geom_segment(aes(x = subject_start, xend = subject_end,
                     y = id2, yend = id2,
                     col = subject), lwd = 1) +
    geom_text(aes(x = (subject_start+subject_end)/2, y = id2, label = sprintf('%0.3f', e_value)),
              size = 2.5, nudge_y = 0.3) +
    geom_segment(data = arm_domains %>% filter(std_name %in% unique(arm_blast$subject[arm_blast$query == q])) %>%
                   mutate(subject = std_name),
                 aes(x = start, xend = end,
                     y = 11, yend = 11),col = 'grey', lwd = 1) +
    scale_x_continuous(breaks = seq(0,5000,500)) +
    scale_y_continuous(breaks = seq(0,100,1),
                       minor_breaks = NULL) +
    labs(x = 'Subject AA Seq') +
    facet_grid(.~subject, space = 'free_x',scale = 'free') +
    theme_bw() +
    theme(axis.title.x = element_text(size = 11),
          axis.text.x = element_text(size = 8),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position = 'None')
  
  qs_plot <- ggarrange(q_plot, s_plot, nrow = 2)
  ggsave(file = sprintf("output/figs/arm_blast_query_%s.jpg", q), qs_plot,
         height = 160, width = 320, units = 'mm',
         bg = 'white',
         dpi = 300)
}

