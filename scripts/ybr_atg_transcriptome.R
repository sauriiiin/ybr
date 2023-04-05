##### YBR196C-A CRISPR ATG MUTANT TRANSCRIPTOME
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 10/14/2022


source(file = 'scripts/initialize.R')

##### COMPILE COUNT DATA
load('~/R/Projects/ybr/output/YBR1_featureCount_raw_sgd_withoverlap.RData')
cnts1 <- as.data.frame(fc$counts[,c(1:3,10:12)])
colnames(cnts1) <- paste('YBR1', c(paste('WT',seq(1,3),sep='_'),paste('ATG',seq(1,3),sep='_')), sep = '_')
head(cnts1)

load('~/R/Projects/ybr/output/YBR3_featureCount_raw_sgd_withoverlap.RData')
cnts2 <- as.data.frame(fc$counts[,c(1:6)])
colnames(cnts2) <- paste('YBR3', c(paste('WT',seq(1,3),sep='_'),paste('ATG',seq(1,3),sep='_')), sep = '_')
head(cnts2)

raw_cnts <- cbind(cnts1, cnts2)
raw_cnts <- raw_cnts[str_detect(rownames(raw_cnts), 'Y') &
                       !str_detect(rownames(raw_cnts), 'tY'),]

##### BATCH CORRECTION
pca_uncorrected_obj = prcomp(raw_cnts)
pca_uncorrected <- as.data.frame(pca_uncorrected_obj[2]$rotation)
pca_uncorrected[,"batch"] <- c(rep('YBR1',6), rep('YBR3',6))
pca_uncorrected[,"samples"] <- rep(c(rep('WT',3), rep('ATG',3)),2)
pca_uncorrected_sum <- summary(pca_uncorrected_obj)

corrected_cnts = ComBat_seq(counts = as.matrix(raw_cnts), batch = c(rep(1,6), rep(2,6)), group = NULL, full_mod = T) #rep(c(rep(1,3), rep(2,3)),2)
pca_corrected_obj = prcomp(corrected_cnts)
pca_corrected <- as.data.frame(pca_corrected_obj[2]$rotation)
pca_corrected[,"batch"] <- c(rep('YBR1',6), rep('YBR3',6))
pca_corrected[,"samples"] <- rep(c(rep('WT',3), rep('ATG',3)),2)
pca_corrected_sum <- summary(pca_corrected_obj)

plot.uncorrected.pca <- pca_uncorrected %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(aes(color=batch, shape=samples), size = 3) +
  scale_color_discrete(name = 'Batch',
                       labels = c('YBR1' = 'One',
                                  'YBR3' = 'Two')) +
  scale_shape_manual(name = 'Strain',
                     values = c('WT' = 19,
                                'ATG' = 17)) +
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
  scale_color_discrete(name = 'Batch',
                       labels = c('YBR1' = 'One',
                                  'YBR3' = 'Two')) +
  scale_shape_manual(name = 'Strain',
                     values = c('WT' = 19,
                                'ATG' = 17)) +
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
# ggsave(sprintf('%sYBR_ATG_SAMPLES_PCA.jpg', fig.path),
#        plot.samples.pca,
#        bg = 'white',
#        width = two.c, height = one.c,
#        units = 'mm',
#        dpi = 300)

##### COUNT CORRELATION
head(raw_cnts)
raw_means <- data.frame(cbind(YBR1_WT = rowSums(raw_cnts[,c(1:3)])/3,
      YBR1_ATG = rowSums(raw_cnts[,c(4:6)])/3,
      YBR3_WT = rowSums(raw_cnts[,c(7:9)])/3,
      YBR3_ATG = rowSums(raw_cnts[,c(10:12)])/3))

corrected_means <- data.frame(cbind(YBR1_WT = rowSums(corrected_cnts[,c(1:3)])/3,
                                    YBR1_ATG = rowSums(corrected_cnts[,c(4:6)])/3,
                                    YBR3_WT = rowSums(corrected_cnts[,c(7:9)])/3,
                                    YBR3_ATG = rowSums(corrected_cnts[,c(10:12)])/3))

plot.b1.rwtvsratg <- raw_means %>%
  ggplot(aes(x = YBR1_WT, y = YBR1_ATG)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman', label.y = 1.8e+05, size = 3) +
  labs(title = 'Raw Batch 1 WT vs ATG',
       x = 'Mean WT Sample Counts',
       y = 'Mean ATG Sample Counts') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt)) +
  coord_cartesian(xlim = c(0,2e+05),
                  ylim = c(0,2e+05))

plot.b2.rwtvsratg <- raw_means %>%
  ggplot(aes(x = YBR3_WT, y = YBR3_ATG)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman', label.y = 1.8e+05, size = 3) +
  labs(title = 'Raw Batch 2 WT vs ATG',
       x = 'Mean WT Sample Counts',
       y = 'Mean ATG Sample Counts') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt)) +
  coord_cartesian(xlim = c(0,2e+05),
                  ylim = c(0,2e+05))

plot.b12.rwtvsrwt <- raw_means %>%
  ggplot(aes(x = YBR1_WT, y = YBR3_WT)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman', label.y = 1.8e+05, size = 3) +
  labs(title = 'Raw Batch 1 WT vs Batch 2 WT',
       x = 'Mean Batch 1 WT Sample Counts',
       y = 'Mean Batch 2 WT Sample Counts') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt)) +
  coord_cartesian(xlim = c(0,2e+05),
                  ylim = c(0,2e+05))

plot.b12.ratgvsratg <- raw_means %>%
  ggplot(aes(x = YBR1_ATG, y = YBR3_ATG)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman', label.y = 1.8e+05, size = 3) +
  labs(title = 'Raw Batch 1 ATG vs Batch 2 ATG',
       x = 'Mean Batch 1 ATG Sample Counts',
       y = 'Mean Batch 2 ATG Sample Counts') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt)) +
  coord_cartesian(xlim = c(0,2e+05),
                  ylim = c(0,2e+05))

plot.b1.rwtvscwt <- corrected_means %>%
  ggplot(aes(x = raw_means$YBR1_WT, y = YBR1_WT)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman', label.y = 1.8e+05, size = 3) +
  labs(title = 'Batch 1 Raw vs Corrected WT',
       x = 'Mean Raw Sample Counts',
       y = 'Mean Corrected Sample Counts') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt)) +
  coord_cartesian(xlim = c(0,2e+05),
                  ylim = c(0,2e+05))

plot.b1.ratgvscatg <- corrected_means %>%
  ggplot(aes(x = raw_means$YBR1_ATG, y = YBR1_ATG)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman', label.y = 1.8e+05, size = 3) +
  labs(title = 'Batch 1 Raw vs Corrected ATG',
       x = 'Mean Raw Sample Counts',
       y = 'Mean Corrected Sample Counts') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt)) +
  coord_cartesian(xlim = c(0,2e+05),
                  ylim = c(0,2e+05))

plot.b2.rwtvscwt <- corrected_means %>%
  ggplot(aes(x = raw_means$YBR3_WT, y = YBR3_WT)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman', label.y = 1.8e+05, size = 3) +
  labs(title = 'Batch 2 Raw vs Corrected WT',
       x = 'Mean Raw Sample Counts',
       y = 'Mean Corrected Sample Counts') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt)) +
  coord_cartesian(xlim = c(0,2e+05),
                  ylim = c(0,2e+05))

plot.b2.ratgvscatg <- corrected_means %>%
  ggplot(aes(x = raw_means$YBR3_ATG, y = YBR3_ATG)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman', label.y = 1.8e+05, size = 3) +
  labs(title = 'Batch 2 Raw vs Corrected ATG',
       x = 'Mean Raw Sample Counts',
       y = 'Mean Corrected Sample Counts') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt)) +
  coord_cartesian(xlim = c(0,2e+05),
                  ylim = c(0,2e+05))

plot.b1.cwtvscatg <- corrected_means %>%
  ggplot(aes(x = YBR1_WT, y = YBR1_ATG)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman', label.y = 1.8e+05, size = 3) +
  labs(title = 'Corrected Batch 1 WT vs ATG',
       x = 'Mean WT Sample Counts',
       y = 'Mean ATG Sample Counts') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt)) +
  coord_cartesian(xlim = c(0,2e+05),
                  ylim = c(0,2e+05))

plot.b2.cwtvscatg <- corrected_means %>%
  ggplot(aes(x = YBR3_WT, y = YBR3_ATG)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman', label.y = 1.8e+05, size = 3) +
  labs(title = 'Corrected Batch 2 WT vs ATG',
       x = 'Mean WT Sample Counts',
       y = 'Mean ATG Sample Counts') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt)) +
  coord_cartesian(xlim = c(0,2e+05),
                  ylim = c(0,2e+05))

plot.b12.cwtvscwt <- corrected_means %>%
  ggplot(aes(x = YBR1_WT, y = YBR3_WT)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman', label.y = 1.8e+05, size = 3) +
  labs(title = 'Corrected Batch 1 WT vs Batch 2 WT',
       x = 'Mean Batch 1 WT Sample Counts',
       y = 'Mean Batch 2 WT Sample Counts') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt)) +
  coord_cartesian(xlim = c(0,2e+05),
                  ylim = c(0,2e+05))

plot.b12.catgvscatg <- corrected_means %>%
  ggplot(aes(x = YBR1_ATG, y = YBR3_ATG)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman', label.y = 1.8e+05, size = 3) +
  labs(title = 'Corrected Batch 1 ATG vs Batch 2 ATG',
       x = 'Mean Batch 1 ATG Sample Counts',
       y = 'Mean Batch 2 ATG Sample Counts') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt)) +
  coord_cartesian(xlim = c(0,2e+05),
                  ylim = c(0,2e+05))


plot.cnt.corr <- plot_grid(plot.b1.rwtvsratg, plot.b2.rwtvsratg, plot.b12.rwtvsrwt, plot.b12.ratgvsratg,
          plot.b1.rwtvscwt, plot.b1.ratgvscatg, plot.b2.rwtvscwt, plot.b2.ratgvscatg,
          plot.b1.cwtvscatg, plot.b2.cwtvscatg, plot.b12.cwtvscwt, plot.b12.catgvscatg,
          nrow = 3, ncol = 4)
plot.cnt.corr 
# ggsave(sprintf('%sYBR_ATG_SAMPLES_COUNTS_CORR.jpg', fig.path),
#        plot.cnt.corr,
#        bg = 'white',
#        width = two.c*2, height = two.c*3/2,
#        units = 'mm',
#        dpi = 300)

##### COMPARISON POST and PRE BATCH CORRECTION
head(raw_cnts)
ybr1_wt_samples = paste(paste('YBR1','WT',sep='_'), seq(1,3), sep = '_')
ybr1_atg_samples = paste(paste('YBR1','ATG',sep='_'), seq(1,3), sep = '_')
ybr3_wt_samples = paste(paste('YBR3','WT',sep='_'), seq(1,3), sep = '_')
ybr3_atg_samples = paste(paste('YBR3','ATG',sep='_'), seq(1,3), sep = '_')
ybr_wt_samples = c(ybr1_wt_samples, ybr3_wt_samples)
ybr_atg_samples = c(ybr1_atg_samples, ybr3_atg_samples)

#run the five comparisons through edgeR using the *uncorrected data*
ybr1_wt_vs_ybr3_wt_uncorrected = run_edgeR(data=raw_cnts, group_a_name="B1_WT", group_a_samples=ybr1_wt_samples, group_b_name="B2_WT", group_b_samples=ybr3_wt_samples)
ybr1_atg_vs_ybr3_atg_uncorrected = run_edgeR(data=raw_cnts, group_a_name="B1_ATG", group_a_samples=ybr1_atg_samples, group_b_name="B2_ATG", group_b_samples=ybr3_atg_samples)
ybr1_wt_vs_ybr1_atg_uncorrected = run_edgeR(data=raw_cnts, group_a_name="B1_WT", group_a_samples=ybr1_wt_samples, group_b_name="B1_ATG", group_b_samples=ybr1_atg_samples)
ybr3_wt_vs_ybr3_atg_uncorrected = run_edgeR(data=raw_cnts, group_a_name="B2_WT", group_a_samples=ybr3_wt_samples, group_b_name="B2_ATG", group_b_samples=ybr3_atg_samples)
ybr_wt_vs_ybr_atg_uncorrected = run_edgeR(data=raw_cnts, group_a_name="Co_WT", group_a_samples=ybr_wt_samples, group_b_name="Co_ATG", group_b_samples=ybr_atg_samples)

#run the same five comparisons through edgeR using the *batch corrected data*
ybr1_wt_vs_ybr3_wt_corrected = run_edgeR(data=corrected_cnts, group_a_name="B1_WT", group_a_samples=ybr1_wt_samples, group_b_name="B2_WT", group_b_samples=ybr3_wt_samples)
ybr1_atg_vs_ybr3_atg_corrected = run_edgeR(data=corrected_cnts, group_a_name="B1_ATG", group_a_samples=ybr1_atg_samples, group_b_name="B2_ATG", group_b_samples=ybr3_atg_samples)
ybr1_wt_vs_ybr1_atg_corrected = run_edgeR(data=corrected_cnts, group_a_name="B1_WT", group_a_samples=ybr1_wt_samples, group_b_name="B1_ATG", group_b_samples=ybr1_atg_samples)
ybr3_wt_vs_ybr3_atg_corrected = run_edgeR(data=corrected_cnts, group_a_name="B2_WT", group_a_samples=ybr3_wt_samples, group_b_name="B2_ATG", group_b_samples=ybr3_atg_samples)
ybr_wt_vs_ybr_atg_corrected = run_edgeR(data=corrected_cnts, group_a_name="Co_WT", group_a_samples=ybr_wt_samples, group_b_name="Co_ATG", group_b_samples=ybr_atg_samples)

#compare wt and atg samples between two batches
sum(ybr1_wt_vs_ybr3_wt_uncorrected$fdr <= 0.05)
sum(ybr1_wt_vs_ybr3_wt_corrected$fdr <= 0.05)

sum(ybr1_atg_vs_ybr3_atg_uncorrected$fdr <= 0.05)
sum(ybr1_atg_vs_ybr3_atg_corrected$fdr <= 0.05)

#wt vs atg results
wt_vs_atg_dea <- rbind(cbind(rbind(ybr1_wt_vs_ybr1_atg_uncorrected, ybr3_wt_vs_ybr3_atg_uncorrected, ybr_wt_vs_ybr_atg_uncorrected), data = 'raw'),
                       cbind(rbind(ybr1_wt_vs_ybr1_atg_corrected, ybr3_wt_vs_ybr3_atg_corrected, ybr_wt_vs_ybr_atg_corrected), data = 'corrected'))
wt_vs_atg_dea$set <- paste(wt_vs_atg_dea$data, wt_vs_atg_dea$contrast, sep = '_')
row.names(wt_vs_atg_dea) <- NULL

wt_vs_atg_dea$DE[wt_vs_atg_dea$fdr <= 0.05 & wt_vs_atg_dea$lfc > 0] <- 'Up'
wt_vs_atg_dea$DE[wt_vs_atg_dea$fdr <= 0.05 & wt_vs_atg_dea$lfc < 0] <- 'Down'

wt_vs_atg_dea %>%
  filter(!is.na(DE)) %>%
  group_by(data,contrast, DE) %>%
  count()

dbWriteTable(conn, 'RNASEQ_YBRATG_DEA', wt_vs_atg_dea, overwrite = T)

set_comparisons <- data.frame(t(combn(unique(wt_vs_atg_dea$set), 2)))
colnames(set_comparisons) <- c('set1','set2')
head(set_comparisons)

all_de_genes <- wt_vs_atg_dea[!is.na(wt_vs_atg_dea$DE), c('set','orf_name','DE')]
for (i in seq(1,nrow(set_comparisons),1)) {
  temp <- all_de_genes %>%
    filter(set %in% c(set_comparisons$set1[i], set_comparisons$set2[i])) %>%
    group_by(DE, orf_name) %>%
    count() %>% filter(n == 2) %>%
    data.frame()
  if (length(temp) > 0) {
    temp$set <- paste(set_comparisons$set1[i], set_comparisons$set2[i], sep = '_and_')
    all_de_genes <- rbind(all_de_genes, temp[,c('set','orf_name','DE')])
  }
  temp <- NULL
}
temp <- all_de_genes %>%
  filter(set %in% c('raw_B1_WT_vs_B1_ATG_and_raw_B2_WT_vs_B2_ATG', 'corrected_B1_WT_vs_B1_ATG_and_corrected_B2_WT_vs_B2_ATG')) %>%
  group_by(DE, orf_name) %>%
  count() %>% filter(n == 2) %>%
  data.frame()
temp$set <- paste('raw_B1_WT_vs_B1_ATG_and_raw_B2_WT_vs_B2_ATG', 'corrected_B1_WT_vs_B1_ATG_and_corrected_B2_WT_vs_B2_ATG', sep = '_and_')
all_de_genes <- rbind(all_de_genes, temp[,c('set','orf_name','DE')])

unique(all_de_genes$set)
all_de_genes_cnt <- all_de_genes %>%
  group_by(set, DE) %>%
  count() %>% data.frame()
# write.csv(all_de_genes_cnt, file = 'output/ybratg_all_de_genes_counts.csv', row.names = F)
dbWriteTable(conn, 'RNASEQ_YBRATG_DEG_COUNTS', all_de_genes_cnt, overwrite = T)

all_genes <- unique(wt_vs_atg_dea$orf_name)
all_genes <- bitr(all_genes, fromType = "ORF",
                  toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
                  OrgDb = org.Sc.sgd.db)
dbWriteTable(conn, 'RNASEQ_YBRATG_ALL_GENES', all_genes, overwrite = T)

goe <- NULL
kegg <- NULL
for (s in unique(all_de_genes$set)) {
  for (de in unique(all_de_genes$DE[all_de_genes$set == s])) {
    temp.deg <- all_de_genes$orf_name[all_de_genes$set == s & all_de_genes$DE == de]
    temp.deg <- bitr(temp.deg, fromType = "ORF",
                      toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
                      OrgDb = org.Sc.sgd.db)
    temp.goe <- enrichGO(gene          = temp.deg$ENSEMBL,
                         universe      = all_genes$ENSEMBL,
                         OrgDb         = org.Sc.sgd.db,
                         keyType       = "ENSEMBL",
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
    if (dim(temp.goe)[1] == 0) {
      cat(sprintf('There are no GO term enrichment for %s DEGs in %s.\n\n',de,s))
    } else{
      cat(sprintf('Top 5 GO term enrichment of %s DEGs in %s for:\n',de,s))
      cat(sprintf('\t1. Cellular component - %s.\n',
                  paste(head(unique(temp.goe$Description[temp.goe$ONTOLOGY == 'CC']),5),collapse = ', ')))
      cat(sprintf('\t2. Biological process - %s.\n',
                  paste(head(unique(temp.goe$Description[temp.goe$ONTOLOGY == 'BP']),5),collapse = ', ')))
      cat(sprintf('\t3. Molecular function - %s.\n\n',
                  paste(head(unique(temp.goe$Description[temp.goe$ONTOLOGY == 'MF']),5),collapse = ', ')))
      goe <- rbind(goe, data.frame(temp.goe, set = s, DE = de))
    }
    
    temp.kegg <- enrichKEGG(gene         = temp.deg$ENSEMBL,
                            universe     = all_genes$ENSEMBL,
                            organism     = 'sce',
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff  = 0.05)
    if (dim(temp.kegg)[1] == 0) {
      cat(sprintf('There are no KEGG pathway enriched for %s DEGs in %s.\n\n',de,s))
    } else{
      cat(sprintf('Top 10 KEGG pathway enrichment for %s DEGs in %s are:\n%s.\n\n',
                  de,s,
                  paste(head(unique(temp.kegg$Description),10),collapse = ', ')))
      kegg <- rbind(kegg, data.frame(temp.kegg, set = s, DE = de))
    }
  }
}
goe$GeneRatio <- as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,2])
goe$BgRatio <- as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,2])
goe$GO <- paste0(goe$ONTOLOGY, '_', goe$Description)
goe <- goe[order(goe$set, goe$DE, goe$GeneRatio, goe$qvalue, decreasing = T),]

kegg$GeneRatio <- as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,2])
kegg$BgRatio <- as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,2])
kegg <- kegg[order(kegg$set, kegg$DE, kegg$GeneRatio, kegg$qvalue, decreasing = T),]

# write.csv(goe, file = 'output/ybratg_goenrich.csv', row.names = FALSE)
# write.csv(kegg, file = 'output/ybratg_keggenrich.csv', row.names = FALSE)
dbWriteTable(conn, 'RNASEQ_YBRATG_GOE', goe, overwrite = T)
dbWriteTable(conn, 'RNASEQ_YBRATG_KEGG', kegg, overwrite = T)

# all_de_genes <- merge(all_de_genes, wt_vs_atg_dea, by = c('set','orf_name','DE'), all.x = T)
# all_de_genes <- merge(all_de_genes, all_genes, by.x = 'orf_name', by.y = 'ORF', all.x = T)
# all_de_genes <- all_de_genes[order(all_de_genes$set, all_de_genes$DE, abs(all_de_genes$lfc), decreasing = T),]
# write.csv(all_de_genes, file = 'output/ybratg_all_de_genes.csv', row.names = FALSE)
dbWriteTable(conn, 'RNASEQ_YBRATG_DEG', all_de_genes, overwrite = T)

##### LFC CORRELATION BETWEEN BATCHES
plot.batch13.lfc.uncorrected <- merge(ybr1_wt_vs_ybr1_atg_uncorrected, ybr3_wt_vs_ybr3_atg_uncorrected,
      by = 'orf_name') %>%
  ggplot(aes(x = lfc.x, y = lfc.y)) +
  geom_point() +
  geom_point(data = merge(ybr1_wt_vs_ybr1_atg_uncorrected, ybr3_wt_vs_ybr3_atg_uncorrected,
                          by = 'orf_name') %>%
               filter(fdr.x <= 0.05, fdr.y <= 0.05),
             aes(x = lfc.x, y = lfc.y), col = 'red') +
  labs(title = 'LFC with Uncorrected Counts',
       x = 'LFC in Batch 1',
       y = 'LFC in Batch 2') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom') +
  coord_cartesian(xlim = c(-4,4),
                  ylim = c(-4,4))

plot.batch13.lfc.corrected <- merge(ybr1_wt_vs_ybr1_atg_corrected, ybr3_wt_vs_ybr3_atg_corrected,
      by = 'orf_name')  %>%
  ggplot(aes(x = lfc.x, y = lfc.y)) +
  geom_point() +
  geom_point(data = merge(ybr1_wt_vs_ybr1_atg_corrected, ybr3_wt_vs_ybr3_atg_corrected,
                          by = 'orf_name') %>%
               filter(fdr.x <= 0.05, fdr.y <= 0.05),
             aes(x = lfc.x, y = lfc.y), col = 'red') +
  labs(title = 'LFC with Corrected Counts',
       x = 'LFC in Batch 1',
       y = 'LFC in Batch 2') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom') +
  coord_cartesian(xlim = c(-4,4),
                  ylim = c(-4,4))

plot.batch13all.lfc.corrected <- merge(merge(ybr1_wt_vs_ybr1_atg_corrected, ybr3_wt_vs_ybr3_atg_corrected,
            by = 'orf_name'),
      ybr_wt_vs_ybr_atg_corrected, by = 'orf_name') %>%
  ggplot(aes(x = lfc.x, y = lfc.y)) +
  geom_point() +
  geom_point(data = merge(merge(ybr1_wt_vs_ybr1_atg_corrected, ybr3_wt_vs_ybr3_atg_corrected,
                                by = 'orf_name'),
                          ybr_wt_vs_ybr_atg_corrected, by = 'orf_name') %>%
               filter(fdr.x <= 0.05, fdr.y <= 0.05, fdr <= 0.05),
             aes(x = lfc.x, y = lfc.y), col = 'red') +
  labs(title = 'LFC with Corrected Counts\nConstrained by Overall Results',
       x = 'LFC in Batch 1',
       y = 'LFC in Batch 2') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom') +
  coord_cartesian(xlim = c(-4,4),
                  ylim = c(-4,4))


plot.batchall1.lfc.corrected <- merge(ybr_wt_vs_ybr_atg_corrected, ybr1_wt_vs_ybr1_atg_corrected,
      by = 'orf_name')  %>%
  ggplot(aes(x = lfc.x, y = lfc.y)) +
  geom_point() +
  geom_point(data = merge(ybr_wt_vs_ybr_atg_corrected, ybr1_wt_vs_ybr1_atg_corrected,
                          by = 'orf_name') %>%
               filter(fdr.x <= 0.05, fdr.y <= 0.05),
             aes(x = lfc.x, y = lfc.y), col = 'red') +
  labs(title = 'LFC with Corrected Counts',
       x = 'LFC in Both Batch Combined',
       y = 'LFC in Batch 1') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom') +
  coord_cartesian(xlim = c(-4,4),
                  ylim = c(-4,4))


plot.batchall3.lfc.corrected <- merge(ybr_wt_vs_ybr_atg_corrected, ybr3_wt_vs_ybr3_atg_corrected,
      by = 'orf_name')  %>%
  ggplot(aes(x = lfc.x, y = lfc.y)) +
  geom_point() +
  geom_point(data = merge(ybr_wt_vs_ybr_atg_corrected, ybr3_wt_vs_ybr3_atg_corrected,
                          by = 'orf_name') %>%
               filter(fdr.x <= 0.05, fdr.y <= 0.05),
             aes(x = lfc.x, y = lfc.y), col = 'red') +
  labs(title = 'LFC with Corrected Counts',
       x = 'LFC in Both Batch Combined',
       y = 'LFC in Batch 2') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom') +
  coord_cartesian(xlim = c(-4,4),
                  ylim = c(-4,4))

plot.compare.lfc <- plot_grid(plot.batch13.lfc.uncorrected, plot.batch13.lfc.corrected, plot.batch13all.lfc.corrected,
          plot.batchall1.lfc.corrected, plot.batchall3.lfc.corrected,
          nrow = 2, ncol = 3,
          align = 'hv')
plot.compare.lfc 
ggsave(sprintf('%sYBR_ATG_COMPARE_LFC.jpg', fig.path),
       plot.compare.lfc,
       bg = 'white',
       width = two.c, height = two.c*2/3,
       units = 'mm', dpi = 300)


##### GO/KEGG ENRICHMENT OF DEGS
ybr_wt_vs_atg_merge_corrected <- merge(merge(ybr1_wt_vs_ybr1_atg_corrected, ybr3_wt_vs_ybr3_atg_corrected,
                                             by = 'orf_name', suffixes = c('_batch1','_batch2')),
                                       ybr_wt_vs_ybr_atg_corrected, by = 'orf_name') %>% data.frame()

ybr_wt_vs_atg_merge_uncorrected <- merge(merge(ybr1_wt_vs_ybr1_atg_uncorrected, ybr3_wt_vs_ybr3_atg_uncorrected,
                                             by = 'orf_name', suffixes = c('_batch1','_batch2')),
                                       ybr_wt_vs_ybr_atg_uncorrected, by = 'orf_name') %>% data.frame()

ybr_wt_vs_atg_merge <- merge(ybr_wt_vs_atg_merge_uncorrected, ybr_wt_vs_atg_merge_corrected, by = 'orf_name', suffixes = c('_uncorrected','_corrected'))



temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_batch1_corrected <= 0.05,
         lfc_batch1_corrected > 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Up', set = 'batch1_corrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_batch1_corrected <= 0.05,
         lfc_batch1_corrected < 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Down', set = 'batch1_corrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_batch2_corrected <= 0.05,
         lfc_batch2_corrected > 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Up', set = 'batch2_corrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_batch2_corrected <= 0.05,
         lfc_batch2_corrected < 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Down', set = 'batch2_corrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_batch1_uncorrected <= 0.05,
         lfc_batch1_uncorrected > 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Up', set = 'batch1_uncorrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_batch1_uncorrected <= 0.05,
         lfc_batch1_uncorrected < 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Down', set = 'batch1_uncorrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_batch2_uncorrected <= 0.05,
         lfc_batch2_uncorrected > 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Up', set = 'batch2_uncorrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_batch2_uncorrected <= 0.05,
         lfc_batch2_uncorrected < 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Down', set = 'batch2_uncorrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_batch1_corrected <= 0.05, fdr_batch2_corrected <= 0.05, fdr_corrected <= 0.05,
         lfc_batch1_corrected > 0, lfc_batch2_corrected > 0, lfc_corrected > 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Up', set = 'alloverlap_corrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_batch1_corrected <= 0.05, fdr_batch2_corrected <= 0.05, fdr_corrected <= 0.05,
         lfc_batch1_corrected < 0, lfc_batch2_corrected < 0, lfc_corrected < 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Down', set = 'alloverlap_corrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_batch1_corrected <= 0.05, fdr_batch2_corrected <= 0.05,
         lfc_batch1_corrected > 0, lfc_batch2_corrected > 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Up', set = 'batch1&2_corrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_batch1_corrected <= 0.05, fdr_batch2_corrected <= 0.05,
         lfc_batch1_corrected < 0, lfc_batch2_corrected < 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Down', set = 'batch1&2_corrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_corrected <= 0.05, lfc_corrected > 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Up', set = 'combined_corrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_corrected <= 0.05, lfc_corrected < 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Down', set = 'combined_corrected'))


temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_batch1_uncorrected <= 0.05, fdr_batch2_uncorrected <= 0.05, fdr_uncorrected <= 0.05,
         lfc_batch1_uncorrected > 0, lfc_batch2_uncorrected > 0, lfc_uncorrected > 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Up', set = 'alloverlap_uncorrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_batch1_uncorrected <= 0.05, fdr_batch2_uncorrected <= 0.05, fdr_uncorrected <= 0.05,
         lfc_batch1_uncorrected < 0, lfc_batch2_uncorrected < 0, lfc_uncorrected < 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
# all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Down', set = 'alloverlap_uncorrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_batch1_uncorrected <= 0.05, fdr_batch2_uncorrected <= 0.05,
         lfc_batch1_uncorrected > 0, lfc_batch2_uncorrected > 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Up', set = 'batch1&2_uncorrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_batch1_uncorrected <= 0.05, fdr_batch2_uncorrected <= 0.05,
         lfc_batch1_uncorrected < 0, lfc_batch2_uncorrected < 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Down', set = 'batch1&2_uncorrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_uncorrected <= 0.05, lfc_uncorrected > 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Up', set = 'combined_uncorrected'))

temp <- ybr_wt_vs_atg_merge %>%
  filter(fdr_uncorrected <= 0.05, lfc_uncorrected < 0)
temp <- bitr(temp$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)
# all_de_genes <- rbind(all_de_genes, cbind(temp, DE = 'Down', set = 'combined_uncorrected'))





head(ybr_wt_vs_atg_merge)
hello <- ybr_wt_vs_atg_merge %>%
  filter(fdr_batch1_uncorrected <= 0.05, fdr_batch2_uncorrected <= 0.05,
         fdr_batch1_corrected <= 0.05, fdr_batch2_corrected <= 0.05,
         lfc_batch1_uncorrected < 0, lfc_batch2_uncorrected < 0,
         lfc_batch1_corrected < 0, lfc_batch2_corrected < 0) 
paste(hello$orf_name, collapse = '; ')


##### COEXPRESSION COIP
ybr.coip <- read.csv('~/R/Projects/ybr/output/ybr_pd_scores.csv')
colnames(ybr.coip)
ybr.coip <- ybr.coip[,c(-1,-7)]
ybr.coip$coip_rank <- as.numeric(rownames(ybr.coip))
dbWriteTable(conn, 'YBR_COIP_RESULTS', ybr.coip)

ybr.coex <- readRDS('/home/aar75/coexpression/20220607_ybr_coexpression.RDS')


ybr.coexip <- merge(ybr.coex, ybr.coip[,-1],
                    by.x = 'gene',
                    by.y = 'ORF',
                    all = T)
ybr.coexip <- ybr.coexip[order(ybr.coexip$rho, decreasing = T),]
rownames(ybr.coexip) <- NULL
ybr.coexip %>%
  filter(transcript %in% c('chr2_614215', 'chr2_614702', 'chr2_614936'))


ybr.coexip.ex <- merge(ybr_wt_vs_atg_merge_corrected, ybr.coexip, by.x = 'orf_name', by.y = 'gene', all = T)
head(ybr.coexip.ex)


## Coip and ybr atg exp
n_all_genes <- ybr.coexip.ex[str_detect(ybr.coexip.ex$orf_name, 'Y') &
                               !str_detect(ybr.coexip.ex$orf_name, 'tY'),] %>%
  filter(!is.na(lfc)) %>%
  nrow()

n_all_coipd_genes <- ybr.coexip.ex[str_detect(ybr.coexip.ex$orf_name, 'Y') &
                                     !str_detect(ybr.coexip.ex$orf_name, 'tY'),] %>%
  filter(!is.na(lfc), !is.na(coip_rank)) %>%
  nrow()

n_dwn_genes <- ybr.coexip.ex[str_detect(ybr.coexip.ex$orf_name, 'Y') &
                               !str_detect(ybr.coexip.ex$orf_name, 'tY'),] %>%
  filter(!is.na(lfc), lfc < 0, fdr <= 0.05) %>%
  nrow()

n_dwn_coipd_genes <- ybr.coexip.ex[str_detect(ybr.coexip.ex$orf_name, 'Y') &
                                     !str_detect(ybr.coexip.ex$orf_name, 'tY'),] %>%
  filter(!is.na(lfc), !is.na(coip_rank), lfc < 0, fdr <= 0.05) %>%
  nrow()

n_up_genes <- ybr.coexip.ex[str_detect(ybr.coexip.ex$orf_name, 'Y') &
                               !str_detect(ybr.coexip.ex$orf_name, 'tY'),] %>%
  filter(!is.na(lfc), lfc > 0, fdr <= 0.05) %>%
  nrow()

n_up_coipd_genes <- ybr.coexip.ex[str_detect(ybr.coexip.ex$orf_name, 'Y') &
                                     !str_detect(ybr.coexip.ex$orf_name, 'tY'),] %>%
  filter(!is.na(lfc), !is.na(coip_rank), lfc > 0, fdr <= 0.05) %>%
  nrow()


coipdwnDETable <- matrix(c(n_dwn_coipd_genes, n_dwn_genes - n_dwn_coipd_genes,
                        n_all_coipd_genes - n_dwn_coipd_genes, n_all_genes + n_dwn_coipd_genes - (n_dwn_genes + n_all_coipd_genes)),
                      nrow = 2,
                      dimnames = list(CoIP=c("yes","no"),
                                      DWN_REG=c("yes","no")))

coipdwnDETable
fisher.test(coipdwnDETable, alternative = "greater")
# coipd proteins are not more likely to be down regulated than not

coipDETable <- matrix(c(n_dwn_coipd_genes + n_up_coipd_genes, n_dwn_genes + n_up_genes - (n_dwn_coipd_genes + n_up_coipd_genes),
                        n_all_coipd_genes - (n_dwn_coipd_genes + n_up_coipd_genes),
                        n_all_genes - n_all_coipd_genes - (n_dwn_genes + n_up_genes - (n_dwn_coipd_genes + n_up_coipd_genes))),
                      nrow = 2,
                      dimnames = list(CoIP=c("yes","no"),
                                      DEG=c("yes","no")))

coipDETable
fisher.test(coipDETable, alternative = "greater")
# coipd proteins are not more likely to be DEd than not

coipDE2Table <- matrix(c(n_dwn_coipd_genes, n_dwn_genes,
                        n_up_coipd_genes, n_up_genes),
                         nrow = 2,
                         dimnames = list(CoIP=c("yes","no"),
                                         DEG=c("dwn","up")))

coipDE2Table
fisher.test(coipDE2Table, alternative = "greater")
# coipd proteins are more likely to be down regulated than up

## deg coipd genes
write.csv(ybr.coexip.ex[str_detect(ybr.coexip.ex$orf_name, 'Y') &
                !str_detect(ybr.coexip.ex$orf_name, 'tY'),] %>%
  filter(!is.na(lfc), !is.na(coip_rank), fdr <= 0.05),
  file = 'output/ybratg_coip_deg.csv', row.names = FALSE)

## Plotting DE counts vs Rho
ybr.coexip.ex$DE[ybr.coexip.ex$fdr <= 0.05 & ybr.coexip.ex$lfc > 0] <- 'Up'
ybr.coexip.ex$DE[ybr.coexip.ex$fdr <= 0.05 & ybr.coexip.ex$lfc < 0] <- 'Down'
ybr.coexip.ex$DE <- factor(ybr.coexip.ex$DE, levels = c('Up','Down'))

plot.rho.vs.degcount <- ybr.coexip.ex %>%
  filter(fdr <= 0.05) %>%
  ggplot( aes(x=rho, fill=DE, group=DE)) +
  geom_histogram(aes(y=..count..),
                 data = ~ subset(., !DE %in% c("Down")),
                 col = 'black',
                 bins = 28) +
  geom_histogram(aes(y=-..count..),
                 data = ~ subset(., DE %in% c("Down")),
                 col = 'black',
                 bins = 28) +
  geom_vline(xintercept = melt(quantile(ybr.coex$rho, c(0.5,0.95,0.99), na.rm = T))[[1]],
             linetype = 'dashed', lwd = 0.4, col = 'red') +
  scale_y_continuous(breaks = seq(-40,40,10)) +
  scale_fill_manual(name = 'DE',
                    breaks = c('Up','Down'),
                    values = c('Up' = '#FFC107',
                               'Down' = '#303F9F')) +
  labs(x = 'Rho',
       y = 'DEG Counts') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom') +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(-30,40))
plot.rho.vs.degcount
ggsave(sprintf('%sYBR_ATG_RHO_VS_DEGCOUNT.jpg', fig.path),
       plot.rho.vs.degcount,
       bg = 'white',
       width = one.5c, height = one.5c,
       units = 'mm', dpi = 300)

## Plotting coexp and ybr atg exp
ybr_atg_rho_vs_lfc_all <- ybr.coexip.ex %>%
  # filter(fdr <= 0.05) %>%
  ggplot(aes(x = rho, y = lfc)) +
  geom_point() +
  geom_point(data = ybr.coexip.ex %>%
               filter(!is.na(coip_rank)),
             aes(x = rho, y = lfc, col = "CoIP'd with YBR"),
             shape = 19) +
  # geom_vline(xintercept = melt(quantile(ybr.coex$rho, c(0.5,0.95,0.99), na.rm = T))[[1]],
  #            linetype = 'dashed', lwd = 0.4, col = 'red') +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman', size = 3, label.x = 0.2, label.y = 2) +
  labs(title = 'Log2FoldChange Inversely Correlates\nwith Coexpression (All Genes)',
       x = 'Rho',
       y = 'LFC in Both Batch Combined') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.position = 'bottom') +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(-2,2))

ybr_atg_rho_vs_lfc_degs <- ybr.coexip.ex %>%
  filter(fdr <= 0.05) %>%
  ggplot(aes(x = rho, y = lfc)) +
  geom_point() +
  geom_point(data = ybr.coexip.ex %>%
               filter(!is.na(coip_rank), fdr <= 0.05),
             aes(x = rho, y = lfc, col = "CoIP'd with YBR"),
             shape = 19) +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman', size = 3, label.x = 0.2, label.y = 2) +
  labs(title = 'Log2FoldChange Inversely Correlates\nwith Coexpression (DEGs)',
       x = 'Rho',
       y = 'LFC in Both Batch Combined') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.position = 'bottom') +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(-2,2))

plot.rho.vs.lfc <- ggarrange(ybr_atg_rho_vs_lfc_all, ybr_atg_rho_vs_lfc_degs,
                                nrow = 1, common.legend = T, legend = 'bottom')
plot.rho.vs.lfc
ggsave(sprintf('%sYBR_ATG_RHO_VS_LFC.jpg', fig.path),
       plot.rho.vs.lfc,
       bg = 'white',
       width = two.c, height = one.c,
       units = 'mm', dpi = 300)







