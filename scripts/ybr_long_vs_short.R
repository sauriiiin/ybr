##### YBR196C-A LONG READ TRANSCRIPTOME
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 11/11/2022


source(file = 'scripts/initialize.R')

##### COMPILE COUNT DATA
load('~/R/Projects/ybr/output/YBR_longread_featureCount_sgd_withoverlap.RData')
cnts1 <- as.data.frame(fc$counts)
colnames(cnts1) <- paste('YBRLR', c(paste('DEL',1,sep='_'),paste('WT',1,sep='_'),paste('STOP',1,sep='_'),paste('SYN',1,sep='_')), sep = '_')
cnts1 <- cnts1[,paste('YBRLR', c(paste('WT',1,sep='_'),paste('DEL',1,sep='_'),paste('SYN',1,sep='_'),paste('STOP',1,sep='_')), sep = '_')]
head(cnts1)

load('~/R/Projects/ybr/output/YBR3_featureCount_raw_sgd_withoverlap.RData')
cnts2 <- as.data.frame(fc$counts[,c(7:12)])
colnames(cnts2) <- paste('YBRSR', c(paste('WT',seq(1,3),sep='_'),paste('DEL',seq(1,3),sep='_')), sep = '_')
head(cnts2)

raw_cnts <- cbind(cnts1, cnts2)
raw_cnts <- raw_cnts[str_detect(rownames(raw_cnts), 'Y') &
                       !str_detect(rownames(raw_cnts), 'tY'),]


##### BATCH CORRECTION
pca_uncorrected_obj = prcomp(raw_cnts)
pca_uncorrected <- as.data.frame(pca_uncorrected_obj[2]$rotation)
pca_uncorrected[,"batch"] <- c(rep('YBRLR',4), rep('YBRSR',6))
pca_uncorrected[,"samples"] <- c(c('WT','DEL','SYN','STOP'),c(rep('WT',3), rep('DEL',3)))
pca_uncorrected_sum <- summary(pca_uncorrected_obj)

corrected_cnts = ComBat_seq(counts = as.matrix(raw_cnts), batch = c(rep(1,4), rep(2,6)), group = NULL, full_mod = T)
pca_corrected_obj = prcomp(corrected_cnts)
pca_corrected <- as.data.frame(pca_corrected_obj[2]$rotation)
pca_corrected[,"batch"] <- c(rep('YBRLR',4), rep('YBRSR',6))
pca_corrected[,"samples"] <- c(c('WT','DEL','SYN','STOP'),c(rep('WT',3), rep('DEL',3)))
pca_corrected_sum <- summary(pca_corrected_obj)

plot.uncorrected.pca <- pca_uncorrected %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(aes(color=batch, shape=samples), size = 3) +
  # scale_color_discrete(name = 'Batch',
  #                      labels = c('YBR1' = 'One',
  #                                 'YBR3' = 'Two')) +
  # scale_shape_manual(name = 'Strain',
  #                    values = c('WT' = 19,
  #                               'DEL' = 17)) +
  labs(title = 'Uncorrected Counts',
       x = sprintf('PC1 (%0.2f%%)', pca_uncorrected_sum$importance[2,1] * 100),
       y = sprintf('PC2 (%0.2f%%)', pca_uncorrected_sum$importance[2,2] * 100)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom')

plot.corrected.pca <- pca_corrected %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(aes(color=batch, shape=samples), size = 3) +
  # scale_color_discrete(name = 'Batch',
  #                      labels = c('YBR1' = 'One',
  #                                 'YBR3' = 'Two')) +
  # scale_shape_manual(name = 'Strain',
  #                    values = c('WT' = 19,
  #                               'DEL' = 17)) +
  labs(title = 'Corrected Counts',
       x = sprintf('PC1 (%0.2f%%)', pca_corrected_sum$importance[2,1] * 100),
       y = sprintf('PC2 (%0.2f%%)', pca_corrected_sum$importance[2,2] * 100)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom')

plot.samples.pca <- ggarrange(plot.uncorrected.pca, plot.corrected.pca,
                              nrow = 1, align = 'hv',
                              common.legend = T, legend = 'bottom')
plot.samples.pca


##### COUNT CORRELATION
corrected_means <- data.frame(cbind(YBRLR_WT = corrected_cnts[,'YBRLR_WT_1'],
                                    YBRLR_DEL = corrected_cnts[,'YBRLR_DEL_1'],
                                    YBRLR_SYN = corrected_cnts[,'YBRLR_SYN_1'],
                                    YBRLR_STOP = corrected_cnts[,'YBRLR_STOP_1'],
                                    YBRSR_WT = rowSums(corrected_cnts[,c(5:7)])/3,
                                    YBRSR_DEL = rowSums(corrected_cnts[,c(8:10)])/3))

corrected_means2 <- corrected_means %>%
  mutate(orf_name = row.names(.)) %>%
  melt(id.vars = 'orf_name', variable.name = 'set', value.name = 'count')

set_combn <- t(combn(unique(corrected_means2$set), 2))
corrected_means_comp <- NULL
for (i in seq(1,nrow(set_combn))) {
  temp <- merge(corrected_means2 %>% filter(set == set_combn[i,1]),
        corrected_means2 %>% filter(set == set_combn[i,2]),
        by = 'orf_name', suffixes = c('_1','_2'))
  corrected_means_comp <- rbind(corrected_means_comp, temp)
}

corrected_means_comp %>%
  ggplot(aes(x = count_1, y = count_2)) +
  geom_point() +
  stat_cor(method = 'spearman') +
  facet_wrap(set_1 ~ set_2, scales = 'free')


##### COMPARISON POST BATCH CORRECTION
ybrlr_wt_samples = paste(paste('YBRLR','WT',sep='_'), 1, sep = '_')
ybrlr_del_samples = paste(paste('YBRLR','DEL',sep='_'), 1, sep = '_')
ybrlr_syn_samples = paste(paste('YBRLR','SYN',sep='_'), 1, sep = '_')
ybrlr_stop_samples = paste(paste('YBRLR','STOP',sep='_'), 1, sep = '_')
ybrsr_wt_samples = paste(paste('YBRSR','WT',sep='_'), seq(1,3), sep = '_')
ybrsr_del_samples = paste(paste('YBRSR','DEL',sep='_'), seq(1,3), sep = '_')

ybr_wt_samples = c(ybrlr_wt_samples, ybrsr_wt_samples)
ybr_del_samples = c(ybrlr_del_samples, ybrsr_del_samples)

ybrlr_wt_vs_ybrsr_wt_corrected = run_edgeR(data=corrected_cnts,
                                           group_a_name="SR_WT", group_a_samples=ybrsr_wt_samples,
                                           group_b_name="LR_WT", group_b_samples=ybrlr_wt_samples)
ybrlr_del_vs_ybrsr_del_corrected = run_edgeR(data=corrected_cnts,
                                           group_a_name="SR_DEL", group_a_samples=ybrsr_del_samples,
                                           group_b_name="LR_DEL", group_b_samples=ybrlr_del_samples)
ybrlr_del_vs_ybrlr_wt_corrected = run_edgeR(data=corrected_cnts,
                                             group_a_name="LR_WT", group_a_samples=ybrlr_wt_samples,
                                             group_b_name="LR_DEL", group_b_samples=ybrlr_del_samples)
ybrsr_del_vs_ybrsr_wt_corrected = run_edgeR(data=corrected_cnts,
                                            group_a_name="SR_WT", group_a_samples=ybrsr_wt_samples,
                                            group_b_name="SR_DEL", group_b_samples=ybrsr_del_samples)
ybr_del_vs_ybr_wt_corrected = run_edgeR(data=corrected_cnts,
                                            group_a_name="ALL_WT", group_a_samples=ybr_wt_samples,
                                            group_b_name="ALL_DEL", group_b_samples=ybr_del_samples)
ybrlr_syn_vs_ybrlr_wt_corrected = run_edgeR(data=corrected_cnts,
                                            group_a_name="LR_WT", group_a_samples=ybrlr_wt_samples,
                                            group_b_name="LR_SYN", group_b_samples=ybrlr_syn_samples)
ybrlr_stop_vs_ybrlr_syn_corrected = run_edgeR(data=corrected_cnts,
                                            group_a_name="LR_SYN", group_a_samples=ybrlr_syn_samples,
                                            group_b_name="LR_STOP", group_b_samples=ybrlr_stop_samples)

ybrlr_del_vs_ybrlr_wt_raw = run_edgeR(data=raw_cnts,
                                            group_a_name="LR_WT", group_a_samples=ybrlr_wt_samples,
                                            group_b_name="LR_DEL", group_b_samples=ybrlr_del_samples)
ybrlr_syn_vs_ybrlr_wt_raw = run_edgeR(data=raw_cnts,
                                            group_a_name="LR_WT", group_a_samples=ybrlr_wt_samples,
                                            group_b_name="LR_SYN", group_b_samples=ybrlr_syn_samples)
ybrlr_stop_vs_ybrlr_syn_raw = run_edgeR(data=raw_cnts,
                                              group_a_name="LR_SYN", group_a_samples=ybrlr_syn_samples,
                                              group_b_name="LR_STOP", group_b_samples=ybrlr_stop_samples)

ybr_del_vs_wt_dea <- data.frame(rbind(cbind(ybrlr_wt_vs_ybrsr_wt_corrected, data = 'corrected'),
                           cbind(ybrlr_del_vs_ybrsr_del_corrected, data = 'corrected'),
                           cbind(ybrlr_del_vs_ybrlr_wt_corrected, data = 'corrected'),
                           cbind(ybrsr_del_vs_ybrsr_wt_corrected, data = 'corrected'),
                           cbind(ybr_del_vs_ybr_wt_corrected, data = 'corrected'),
                           cbind(ybrlr_syn_vs_ybrlr_wt_corrected, data = 'corrected'),
                           cbind(ybrlr_stop_vs_ybrlr_syn_corrected, data = 'corrected'),
                           cbind(ybrlr_del_vs_ybrlr_wt_raw, data = 'raw'),
                           cbind(ybrlr_syn_vs_ybrlr_wt_raw, data = 'raw'),
                           cbind(ybrlr_stop_vs_ybrlr_syn_raw, data = 'raw')))
row.names(ybr_del_vs_wt_dea) <- NULL

ybr_del_vs_wt_dea$DE[ybr_del_vs_wt_dea$fdr <= 0.05 & ybr_del_vs_wt_dea$lfc > 0] <- 'Up'
ybr_del_vs_wt_dea$DE[ybr_del_vs_wt_dea$fdr <= 0.05 & ybr_del_vs_wt_dea$lfc < 0] <- 'Down'


ybr_del_vs_wt_dea %>%
  filter(!is.na(DE)) %>%
  group_by(data, contrast, DE) %>%
  dplyr::count()

ybr_del_vs_wt_dea %>%
  filter(contrast == 'LR_WT_vs_LR_DEL', data == 'raw', !is.na(DE))


dbWriteTable(conn, 'RNASEQ_YBRSRLR_DEA', ybr_del_vs_wt_dea, overwrite = T)
