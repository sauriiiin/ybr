
source(file = 'scripts/initialize.R')


ybratg.dea <- dbGetQuery(conn, 'select * from RNASEQ_YBRATG_DEA')
head(ybratg.dea)
unique(ybratg.dea$set)

plot.ybratg.dea.raw <- ybratg.dea %>%
  filter(data == 'raw') %>%
  ggplot(aes(x = lfc, y = -log10p, col = DE)) +
  geom_point(size = 0.5) +
  labs(x = 'logFoldChange') +
  facet_wrap(.~set, labeller = labeller(set = c('raw_B1_WT_vs_B1_ATG' = 'Raw Batch1\nWT vs ATG',
                                                'raw_B2_WT_vs_B2_ATG' = 'Raw Batch2\nWT vs ATG',
                                                'raw_Co_WT_vs_Co_ATG' = 'Raw Combined\nWT vs ATG',
                                                'corrected_B1_WT_vs_B1_ATG' = 'Corrected Batch1\nWT vs ATG',
                                                'corrected_B2_WT_vs_B2_ATG' = 'Corrected Batch2\nWT vs ATG',
                                                'corrected_Co_WT_vs_Co_ATG' = 'Correcred Combined\nWT vs ATG'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom') +
  coord_cartesian(xlim = c(-4,4),
                  ylim = c(0,150))
ggsave(sprintf('%sYBR_ATG_RAW_VOLCANO.jpg', fig.path),
       plot.ybratg.dea.raw,
       bg = 'white',
       width = two.c, height = one.c,
       units = 'mm',
       dpi = 300)


plot.ybratg.dea.cor <- ybratg.dea %>%
  filter(data == 'corrected') %>%
  ggplot(aes(x = lfc, y = -log10p, col = DE)) +
  geom_point(size = 0.5) +
  labs(x = 'logFoldChange') +
  facet_wrap(.~set, labeller = labeller(set = c('raw_B1_WT_vs_B1_ATG' = 'Raw Batch1\nWT vs ATG',
                                                'raw_B2_WT_vs_B2_ATG' = 'Raw Batch2\nWT vs ATG',
                                                'raw_Co_WT_vs_Co_ATG' = 'Raw Combined\nWT vs ATG',
                                                'corrected_B1_WT_vs_B1_ATG' = 'Corrected Batch1\nWT vs ATG',
                                                'corrected_B2_WT_vs_B2_ATG' = 'Corrected Batch2\nWT vs ATG',
                                                'corrected_Co_WT_vs_Co_ATG' = 'Correcred Combined\nWT vs ATG'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom') +
  coord_cartesian(xlim = c(-4,4),
                  ylim = c(0,150))
ggsave(sprintf('%sYBR_ATG_CORRECTED_VOLCANO.jpg', fig.path),
       plot.ybratg.dea.cor,
       bg = 'white',
       width = two.c, height = one.c,
       units = 'mm',
       dpi = 300)


#####
head(ybratg.dea)

temp <- merge(ybratg.dea %>% filter(set == 'raw_B1_WT_vs_B1_ATG'),
              ybratg.dea %>% filter(set == 'raw_B2_WT_vs_B2_ATG'),
              by = c('data','orf_name'), suffixes = c('_batch1','_batch2'))
plot.lfc.b1b2 <- temp %>%
  ggplot(aes(x = lfc_batch1, y = lfc_batch2)) +
  geom_point() +
  geom_point(data = temp %>% filter(!is.na(DE_batch1), !is.na(DE_batch1), DE_batch1 == DE_batch2),
             aes(x = lfc_batch1, y = lfc_batch2), col = 'red') +
  labs(x = 'logFoldChange in Batch 1',
       y = 'logFoldChange in Batch 2') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom') +
  coord_cartesian(xlim = c(-4,4),
                  ylim = c(-4,4))
ggsave(sprintf('%sYBR_ATG_RAW_LFC_COMPARE_B1B2.jpg', fig.path),
       plot.lfc.b1b2,
       bg = 'white',
       width = one.c, height = one.c,
       units = 'mm',
       dpi = 300)

temp <- merge(ybratg.dea %>% filter(set == 'corrected_B1_WT_vs_B1_ATG'),
              ybratg.dea %>% filter(set == 'corrected_B2_WT_vs_B2_ATG'),
              by = c('data','orf_name'), suffixes = c('_batch1','_batch2'))
plot.cor.lfc.b1b2 <- temp %>%
  ggplot(aes(x = lfc_batch1, y = lfc_batch2)) +
  geom_point() +
  geom_point(data = temp %>% filter(!is.na(DE_batch1), !is.na(DE_batch1), DE_batch1 == DE_batch2),
             aes(x = lfc_batch1, y = lfc_batch2), col = 'red') +
  labs(x = 'logFoldChange in Batch 1',
       y = 'logFoldChange in Batch 2') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom') +
  coord_cartesian(xlim = c(-4,4),
                  ylim = c(-4,4))
ggsave(sprintf('%sYBR_ATG_COR_LFC_COMPARE_B1B2.jpg', fig.path),
       plot.cor.lfc.b1b2,
       bg = 'white',
       width = one.c, height = one.c,
       units = 'mm',
       dpi = 300)
