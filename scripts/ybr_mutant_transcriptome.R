##### YBR196C-A MUTANT TRANSCRIPTOME
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 10/27/2022


source(file = 'scripts/initialize.R')

# ##### GATHER ORF DATA
# orfids <- dbGetQuery(conn, 'select b.orf_id, a.gene_systematic_name, a.orf_class, a.is_transient, a.is_candidate, a.informative
#                             from
#                             aaron.scer_orfs_translation_info a, CARVUNIS_YEAST.orf_ids_crossmap b
#                             where a.orf_id = b.aaron_orfid')
# orfids <- orfids[-1,]
# orfids <- orfids %>%
#   group_by(orf_id, gene_systematic_name, orf_class, is_transient, is_candidate, informative) %>%
#   count() %>%
#   data.frame()
# orfids <- orfids[,-7]
# orfids$transient[orfids$is_candidate == 1 & orfids$informative] <- 'No'
# orfids$transient[orfids$is_candidate + orfids$is_transient  == 2] <- 'Yes'
# orfids$transient[is.na(orfids$transient)] <- 'NEI'
# 
# 
# sum(orfids$is_candidate == 1 & orfids$is_transient == 0)
# 
# orfids %>%
#   group_by(transient) %>%
#   count

##### COMPILE COUNT DATA
load('~/R/Projects/ybr/output/YBR1_featureCount_raw_sgd_withoverlap.RData')
cnts1 <- as.data.frame(fc$counts[,c(1:3,7:12)])
colnames(cnts1) <- paste('YBR1', c(paste('WT',seq(1,3),sep='_'),paste('DEL',seq(1,3),sep='_'),paste('ATG',seq(1,3),sep='_')), sep = '_')
head(cnts1)

load('~/R/Projects/ybr/output/YBR3_featureCount_raw_sgd_withoverlap.RData')
cnts2 <- as.data.frame(fc$counts[,c(1:6)])
# colnames(cnts2) <- paste('YBR3', c(paste('WT',seq(1,3),sep='_'),paste('ATG',seq(1,3),sep='_'),
#                                    paste('WT2',seq(1,3),sep='_'),paste('DEL',seq(1,3),sep='_')), sep = '_')
colnames(cnts2) <- paste('YBR3', c(paste('WT',seq(1,3),sep='_'),paste('ATG',seq(1,3),sep='_')), sep = '_')
head(cnts2)

raw_cnts <- cbind(cnts1, cnts2)
raw_cnts <- raw_cnts[str_detect(rownames(raw_cnts), 'Y') &
                       !str_detect(rownames(raw_cnts), 'tY'),]

##### BATCH CORRECTION
pca_uncorrected_obj = prcomp(raw_cnts)
pca_uncorrected <- as.data.frame(pca_uncorrected_obj[2]$rotation)
# pca_uncorrected[,"batch"] <- c(rep('YBR1',9), rep('YBR3',12))
# pca_uncorrected[,"bg"] <- c(rep('BY4742',15), rep('BY4741',6))
pca_uncorrected[,"batch"] <- c(rep('YBR1',9), rep('YBR3',6))
pca_uncorrected[,"bg"] <- c(rep('BY4742',15))
pca_uncorrected[,"samples"] <- str_split(row.names(pca_uncorrected), '_', simplify = T)[,2]
pca_uncorrected_sum <- summary(pca_uncorrected_obj)

# corrected_cnts = ComBat_seq(counts = as.matrix(raw_cnts), batch = c(rep(1,9), rep(2,12)), group = c(rep(1,15), rep(2,6)), full_mod = T)
corrected_cnts = ComBat_seq(counts = as.matrix(raw_cnts), batch = c(rep(1,9), rep(2,6)), group = NULL, full_mod = T)
pca_corrected_obj = prcomp(corrected_cnts)
pca_corrected <- as.data.frame(pca_corrected_obj[2]$rotation)
# pca_corrected[,"batch"] <- c(rep('YBR1',9), rep('YBR3',12))
# pca_corrected[,"bg"] <- c(rep('BY4742',15), rep('BY4741',6))
pca_corrected[,"batch"] <- c(rep('YBR1',9), rep('YBR3',6))
pca_corrected[,"bg"] <- c(rep('BY4742',15))
pca_corrected[,"samples"] <- str_split(row.names(pca_corrected), '_', simplify = T)[,2]
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
  facet_wrap(.~bg) +
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
  facet_wrap(.~bg) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom')

plot.samples.pca <- ggarrange(plot.uncorrected.pca, plot.corrected.pca,
                              nrow = 2, align = 'hv',
                              common.legend = T, legend = 'bottom')
plot.samples.pca

##### MEAN COUNTS
raw_cnts_sum <- raw_cnts %>%
  mutate(orf_name = row.names(.)) %>%
  melt(id.vars = 'orf_name', variable.name = 'sid', value.name = 'count')
temp <- cbind(unique(as.character(raw_cnts_sum$sid)), str_split(unique(raw_cnts_sum$sid), '_', simplify = T)) %>% data.frame()
colnames(temp) <- c('sid','batch','sample','replicate')
raw_cnts_sum <- merge(temp, raw_cnts_sum, by = 'sid', all = T)
raw_cnts_sum <- raw_cnts_sum %>%
  group_by(batch, sample, orf_name) %>%
  summarize(mean_cnt = mean(count, na.rm = T),
            median_cnt = median(count, na.rm = T),
            std_cnt = sd(count, na.rm = T),
            .groups = 'keep') %>%
  mutate(cv_cnt = std_cnt/mean_cnt) %>%
  data.frame()
head(raw_cnts_sum)

corrected_cnts_sum <- corrected_cnts %>%
  data.frame() %>%
  mutate(orf_name = row.names(.)) %>%
  melt(id.vars = 'orf_name', variable.name = 'sid', value.name = 'count')
temp <- cbind(unique(as.character(corrected_cnts_sum$sid)), str_split(unique(corrected_cnts_sum$sid), '_', simplify = T)) %>% data.frame()
colnames(temp) <- c('sid','batch','sample','replicate')
corrected_cnts_sum <- merge(temp, corrected_cnts_sum, by = 'sid', all = T)
corrected_cnts_sum <- corrected_cnts_sum %>%
  group_by(batch, sample, orf_name) %>%
  summarize(mean_cnt = mean(count, na.rm = T),
            median_cnt = median(count, na.rm = T),
            std_cnt = sd(count, na.rm = T),
            .groups = 'keep') %>%
  mutate(cv_cnt = std_cnt/mean_cnt) %>%
  data.frame()
head(corrected_cnts_sum)

corrected_cnts_sum %>%
  filter(batch == 'YBR1', sample == 'WT') %>%
  ggplot(aes(x = median_cnt, y = corrected_cnts_sum$median_cnt[corrected_cnts_sum$batch == 'YBR1' &
                                                         corrected_cnts_sum$sample == 'ATG'])) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman')


##### DEA
ybr1_wt_samples = paste(paste('YBR1','WT',sep='_'), seq(1,3), sep = '_')
ybr1_atg_samples = paste(paste('YBR1','ATG',sep='_'), seq(1,3), sep = '_')
ybr1_del_samples = paste(paste('YBR1','DEL',sep='_'), seq(1,3), sep = '_')

ybr3_wt_samples = paste(paste('YBR3','WT',sep='_'), seq(1,3), sep = '_')
ybr3_atg_samples = paste(paste('YBR3','ATG',sep='_'), seq(1,3), sep = '_')

# ybr3_wt2_samples = paste(paste('YBR3','WT2',sep='_'), seq(1,3), sep = '_')
# ybr3_del_samples = paste(paste('YBR3','DEL',sep='_'), seq(1,3), sep = '_')

ybr_wt_samples = c(ybr1_wt_samples, ybr3_wt_samples)
ybr_atg_samples = c(ybr1_atg_samples, ybr3_atg_samples)


ybr1_wt_vs_ybr3_wt_corrected = run_edgeR(data=corrected_cnts, group_a_name="B1_WT", group_a_samples=ybr1_wt_samples, group_b_name="B2_WT", group_b_samples=ybr3_wt_samples)
ybr1_atg_vs_ybr3_atg_corrected = run_edgeR(data=corrected_cnts, group_a_name="B1_ATG", group_a_samples=ybr1_atg_samples, group_b_name="B2_ATG", group_b_samples=ybr3_atg_samples)
ybr1_wt_vs_ybr1_atg_corrected = run_edgeR(data=corrected_cnts, group_a_name="B1_WT", group_a_samples=ybr1_wt_samples, group_b_name="B1_ATG", group_b_samples=ybr1_atg_samples)
ybr3_wt_vs_ybr3_atg_corrected = run_edgeR(data=corrected_cnts, group_a_name="B2_WT", group_a_samples=ybr3_wt_samples, group_b_name="B2_ATG", group_b_samples=ybr3_atg_samples)
ybr_wt_vs_ybr_atg_corrected = run_edgeR(data=corrected_cnts, group_a_name="Co_WT", group_a_samples=ybr_wt_samples, group_b_name="Co_ATG", group_b_samples=ybr_atg_samples)

ybr1_wt_vs_ybr1_del_corrected = run_edgeR(data=corrected_cnts, group_a_name="B1_WT", group_a_samples=ybr1_wt_samples, group_b_name="B1_DEL", group_b_samples=ybr1_del_samples)
# ybr3_wt2_vs_ybr3_del_corrected = run_edgeR(data=corrected_cnts, group_a_name="B2_WT2", group_a_samples=ybr3_wt2_samples, group_b_name="B2_DEL", group_b_samples=ybr3_del_samples)

wt_vs_mut_dea <- rbind(cbind(rbind(ybr1_wt_vs_ybr1_atg_corrected, ybr3_wt_vs_ybr3_atg_corrected, ybr_wt_vs_ybr_atg_corrected), mutant = 'atg'),
                       cbind(rbind(ybr1_wt_vs_ybr1_del_corrected), mutant = 'del'))
# wt_vs_mut_dea$set <- paste(wt_vs_mut_dea$mutant, wt_vs_mut_dea$contrast, sep = '_')
row.names(wt_vs_mut_dea) <- NULL

wt_vs_mut_dea$DE[wt_vs_mut_dea$fdr <= 0.05 & wt_vs_mut_dea$lfc > 0] <- 'Up'
wt_vs_mut_dea$DE[wt_vs_mut_dea$fdr <= 0.05 & wt_vs_mut_dea$lfc < 0] <- 'Down'

wt_vs_mut_dea %>%
  filter(!is.na(DE)) %>%
  group_by(mutant, contrast, DE) %>%
  count()

dbWriteTable(conn, 'RNASEQ_YBRMUT_DEA', wt_vs_mut_dea, overwrite = T)


listInputUP = list("Batch1 WT vs Batch1 ATG" = ybr1_wt_vs_ybr1_atg_corrected$orf_name[ybr1_wt_vs_ybr1_atg_corrected$fdr <= 0.05 &
                                                                                        ybr1_wt_vs_ybr1_atg_corrected$lfc > 0],
                  "Batch2 WT vs Batch2 ATG" = ybr3_wt_vs_ybr3_atg_corrected$orf_name[ybr3_wt_vs_ybr3_atg_corrected$fdr <= 0.05 &
                                                                                       ybr3_wt_vs_ybr3_atg_corrected$lfc > 0],
                  "Combined WT vs Combined ATG" = ybr_wt_vs_ybr_atg_corrected$orf_name[ybr_wt_vs_ybr_atg_corrected$fdr <= 0.05 &
                                                                                         ybr_wt_vs_ybr_atg_corrected$lfc > 0])

upset(fromList(listInputUP), order.by = "freq", point.size=5, set_size.scale_max = 6000, text.scale = 1)


listInputDN = list("Batch1 WT vs Batch1 ATG" = ybr1_wt_vs_ybr1_atg_corrected$orf_name[ybr1_wt_vs_ybr1_atg_corrected$fdr <= 0.05 &
                                                                                        ybr1_wt_vs_ybr1_atg_corrected$lfc < 0],
                   "Batch2 WT vs Batch2 ATG" = ybr3_wt_vs_ybr3_atg_corrected$orf_name[ybr3_wt_vs_ybr3_atg_corrected$fdr <= 0.05 &
                                                                                        ybr3_wt_vs_ybr3_atg_corrected$lfc < 0],
                   "Combined WT vs Combined ATG" = ybr_wt_vs_ybr_atg_corrected$orf_name[ybr_wt_vs_ybr_atg_corrected$fdr <= 0.05 &
                                                                                          ybr_wt_vs_ybr_atg_corrected$lfc < 0])

upset(fromList(listInputDN), order.by = "freq", point.size=5, set_size.scale_max = 6000, text.scale = 1)



##### LFC BASED GO AND KEGG ENRICHMENT IN COMMON
head(wt_vs_mut_dea)
all_genes <- unique(wt_vs_mut_dea$orf_name)
all_genes <- bitr(all_genes, fromType = "ORF",
                  toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
                  OrgDb = org.Sc.sgd.db)

kegg <- NULL
goe <- NULL
for (c in unique(wt_vs_mut_dea$contrast)) {
  for (lfc in seq(0,3,0.5)) {
    for (de in c(-1,1)) {
      temp.deg <- NULL
      if (de == -1) {
        temp.deg <- wt_vs_mut_dea$orf_name[wt_vs_mut_dea$contrast == c & wt_vs_mut_dea$fdr <= 0.05 &
                                             wt_vs_mut_dea$lfc < lfc*de]
      } else {
        temp.deg <- wt_vs_mut_dea$orf_name[wt_vs_mut_dea$contrast == c & wt_vs_mut_dea$fdr <= 0.05 &
                                             wt_vs_mut_dea$lfc > lfc]
      }
      
      if (length(temp.deg) > 0) {
        cat(sprintf('There are %d DEGs at lfc > %0.2f in %s.\n\n',nrow(temp.deg),lfc*de,c))
        
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
        if (length(temp.goe) > 0) {
          if (dim(temp.goe)[1] == 0) {
            cat(sprintf('There are no GO term enrichment for lfc > %0.2f in %s.\n\n',lfc*de,c))
          } else{
            goe <- rbind(goe, data.frame(temp.goe, contrast = c, lfc = lfc*de))
          }
        }
        
        temp.kegg <- enrichKEGG(gene         = temp.deg$ENSEMBL,
                                universe     = all_genes$ENSEMBL,
                                organism     = 'sce',
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff  = 0.05)
        if (length(temp.kegg) > 0) {
          if (dim(temp.kegg)[1] == 0) {
            cat(sprintf('There are no KEGG pathway enriched for lfc > %0.2f in %s.\n\n',lfc*de,c))
          } else{
            kegg <- rbind(kegg, data.frame(temp.kegg, contrast = c, lfc = lfc*de))
          }
        }
      }
    }
  }
}
goe$GeneRatio <- as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,2])
goe$BgRatio <- as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,2])
goe$GO <- paste0(goe$ONTOLOGY, '_', goe$Description)
goe <- goe[order(goe$contrast, goe$lfc, goe$GeneRatio, goe$qvalue, decreasing = T),]

kegg$GeneRatio <- as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,2])
kegg$BgRatio <- as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,2])
kegg <- kegg[order(kegg$contrast, kegg$lfc, kegg$GeneRatio, kegg$qvalue, decreasing = T),]

head(goe)
goe %>%
  group_by(contrast, lfc) %>%
  count() %>% data.frame()

dbWriteTable(conn, 'RNASEQ_YBRMUT_GOE', goe, overwrite = T)
dbWriteTable(conn, 'RNASEQ_YBRMUT_KEGG', kegg, overwrite = T)


##### COMMON FOR ATG AND DEL
hello <- wt_vs_mut_dea %>%
  filter(contrast %in% c('Co_WT_vs_Co_ATG', 'B1_WT_vs_B1_DEL'), !is.na(DE)) %>%
  group_by(orf_name, DE) %>%
  count() %>%
  filter(n == 2) %>%
  data.frame()
write.csv(hello, file = 'output/ybrmutcommon.csv', row.names = F)


