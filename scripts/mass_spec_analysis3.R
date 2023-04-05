##### MASS SPEC RESULT OVERLAP AND COEXPRESSION DATA
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 01/16/2023

##### INITIALIZE
library(RMariaDB)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(GO.db)
library(dplyr)
library(stringr)
library(reshape2)
source('/home/sbp29/R/Projects/ybr/scripts/run_edgeR.R')
`%notin%` <- Negate(`%in%`)

##### LOAD PRESAVED DATA
load('~/R/Projects/ybr/data/230119_pulldown_rnaseq_coexp_init.RData')
conn <- initialize.sql("saurin_test")

##### LOAD COUNT DATA
pdone.cnts <- read.csv(file = '/home/sbp29/R/Projects/ybr/scripts/pull_down_results/pull_down_one_cnts.csv')
pdtwo.cnts <- read.csv(file = '/home/sbp29/R/Projects/ybr/scripts/pull_down_results/pull_down_two_cnts.csv')
sample_info <- read.csv(file = '/home/sbp29/R/Projects/ybr/input/ms_data/221225_sample_info.csv')
sample_info$Sample.ID.number <- str_replace_all(sample_info$Sample.ID.number, '-', '.')

pddat.cnts <- rbind(melt(pdone.cnts, id.vars = 'Gene', variable.name = 'sample', value.name = 'peptide_cnts'),
                    melt(pdtwo.cnts, id.vars = 'Gene', variable.name = 'sample', value.name = 'peptide_cnts'))
pddat.cnts <- merge(pddat.cnts, sample_info[,c(1,2,8:10)], by.x = 'sample', by.y = 'Sample.ID.number')
pddat.cnts <- pddat.cnts %>%
  mutate(ID = paste('attempt', Attempt, background, ybr_tag, sep = '_'))
pddat.cnts <- pddat.cnts %>%
  group_by(Attempt, background, ybr_tag, replicate, ID, Gene) %>%
  summarise(peptide_cnts = sum(peptide_cnts), .groups = 'keep') %>%
  data.frame()


##### FILTER HA DATA AND MAKE MATRIX
pddat.cnts.ha <- pddat.cnts %>% filter(ID %in% c('attempt_2_WT_YBR_HA','attempt_2_WT_HA'))
col.names <- NULL
pddat.cnts.ha.mat <- data.frame(Gene = unique(pddat.cnts.ha$Gene))
for (id in unique(pddat.cnts.ha$ID)) {
  for (r in unique(pddat.cnts.ha$replicate[pddat.cnts.ha$ID == id])) {
    temp <- pddat.cnts.ha %>% filter(ID == id, replicate == r)
    pddat.cnts.ha.mat <- merge(pddat.cnts.ha.mat, temp[,c('Gene','peptide_cnts')], by = 'Gene', all.x = T)
    col.names <- c(col.names, paste(id, r, sep = '_'))
  }
}
colnames(pddat.cnts.ha.mat) <- c('Gene', col.names)
pddat.cnts.ha.mat[is.na(pddat.cnts.ha.mat)] <- 0

edgeR_in <- pddat.cnts.ha.mat[-1]
row.names(edgeR_in) <- pddat.cnts.ha.mat$Gene
hello <- rowSums(edgeR_in) %>% data.frame() %>% filter(. > 0)
edgeR_in <- edgeR_in[row.names(edgeR_in) %in% row.names(hello),]
  
edgeR_out <- run_edgeR(data=edgeR_in,
                       group_a_name='Control', group_a_samples=col.names[1:3],
                       group_b_name='YBR', group_b_samples=col.names[4:6])
edgeR_out %>% filter(lfc > 0, fdr <= 0.2)
edgeR_out <- edgeR_out[-5]

##### RUN SOME OTHER TEST AS WELL
## NORMALIZE BY LENGTH AND LIBRARY SIZE
pddat.cnts.ha <- merge(pddat.cnts.ha, pddat.cnts.ha %>%
        group_by(ID, replicate) %>%
        summarize(total_peptide_cnts = sum(peptide_cnts, na.rm = T), .groups = 'keep'),
      by = c('ID','replicate')) %>%
  mutate(norm_peptide_cnts = peptide_cnts/total_peptide_cnts)

for (g in unique(edgeR_out$orf_name)) {
  temp <- pddat.cnts.ha %>% filter(Gene == g)
  temp.kruskal <- compare_means(norm_peptide_cnts ~ ID, temp, method = 'kruskal')
  temp.anova <- compare_means(norm_peptide_cnts ~ ID, temp, method = 'anova')
  
  edgeR_out$kruskal[edgeR_out$orf_name == g] <- temp.kruskal$p
  edgeR_out$anova[edgeR_out$orf_name == g] <- temp.anova$p
}
head(edgeR_out)
p.adjust(edgeR_out$anova, 'BH')

##### WRITE THE RESULTS
overlap_gene.descriptions <- all_genes[(all_genes$ORF %in% edgeR_out$orf_name) |
                                         (all_genes$GENENAME %in% edgeR_out$orf_name),]
edgeR_out <- rbind(merge(edgeR_out, overlap_gene.descriptions, by.x = 'orf_name', by.y = 'GENENAME')[c(1:6,11)],
                             merge(edgeR_out, overlap_gene.descriptions, by.x = 'orf_name', by.y = 'ORF')[c(1:6,11)])
edgeR_out <- edgeR_out[order(edgeR_out$lfc, decreasing = T),]
rownames(edgeR_out) <- NULL
write.csv(edgeR_out, file = 'output/pull_down_two_ha_samples_enrichment.csv', row.names = F)


##### OVERLAP WITH MNG SAMPLES
overlap_gene.two_mng <- read.csv(file = 'output/pull_down_two_mNG_samples_overlap.csv')

fdr_lim <- 0.3
edgeR_out %>%
  filter(lfc > 0, fdr <= fdr_lim) %>%
  nrow()

hello <- edgeR_out %>%
  filter(lfc > 0, fdr <= fdr_lim, orf_name %in% overlap_gene.two_mng$Gene[overlap_gene.two_mng$lfc > 0])

paste(hello$orf_name, collapse = ',')

# write.csv(hello, file = 'output/ha_and_mng_fdr05.csv')

##### GSEA OF PULLED DOWN PROTEINS
ybr_ha_ppi <- rbind(merge(edgeR_out[-7] %>% filter(lfc > 0, fdr <= 0.1), overlap_gene.descriptions, by.x = 'orf_name', by.y = 'GENENAME')[c(1:6,10,11)],
      merge(edgeR_out[-7] %>% filter(lfc > 0, fdr <= 0.1), overlap_gene.descriptions, by.x = 'orf_name', by.y = 'ORF')[c(1:6,10,11)])
ybr_ha_ppi <- ybr_ha_ppi[order(ybr_ha_ppi$lfc, decreasing = T),]
head(ybr_ha_ppi)

temp.goe <- enrichGO(gene          = ybr_ha_ppi$ENSEMBL,
                     universe      = all_genes$ENSEMBL,
                     OrgDb         = org.Sc.sgd.db,
                     keyType       = "ENSEMBL",
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)
data.frame(temp.goe)
coip_two_ha_go <- data.frame(temp.goe)
coip_two_ha_go$GeneRatio <- as.numeric(str_split(coip_two_ha_go$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(coip_two_ha_go$GeneRatio,'/',simplify = T)[,2])
coip_two_ha_go$BgRatio <- as.numeric(str_split(coip_two_ha_go$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(coip_two_ha_go$BgRatio,'/',simplify = T)[,2])
write.csv(coip_two_ha_go, file = 'output/figs/ybr_presentation/pull_down_two_ha_edger_GO.csv')

plot.ha.emap <- emapplot(enrichplot::pairwise_termsim(temp.goe), showCategory = 30)
ggsave("output/figs/ybr_presentation/pull_down_two_ha_edger_emap.jpg", plot.ha.emap,
       height = 400, width = 400, units = 'mm',
       bg = 'white',
       dpi = 300)


#####
head(ybr_ha_ppi)
plot.coip_coexp <- coexp.ybr %>%
  filter(!is.na(orf_name)) %>%
  mutate(set = case_when(orf_name %in% ybr_ha_ppi$ENSEMBL ~ "CoIP'd with\nYBR196C-A-HA",
                         TRUE ~ 'All')) %>%
  ggplot(aes(x = set, y = cor)) +
  labs(y = 'Co-expression with YBR196C-A transcript') +
  # geom_violin() +
  # geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75)) +
  stat_compare_means(hjust = -0.2, label.y = 1, size = 3) +
  theme_linedraw() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 8),
        axis.title.y = element_text(size = 10))
ggsave("output/figs/ybr_presentation/pull_down_two_ha_coexp.jpg", plot.coip_coexp,
       height = 100, width = 100, units = 'mm',
       bg = 'white',
       dpi = 300)
