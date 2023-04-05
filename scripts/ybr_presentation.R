
gcn4_targets <- read.csv(file = 'input/GCN4_targets2.csv')

head(dea_results)
unique(dea_results$contrast)

ybr_del_dea <- dea_results[dea_results$contrast == 'YBR3_BY4741_DEL_vs_YBR3_BY4741_WT',]
ybr_del_dea <- ybr_del_dea[order(ybr_del_dea$lfc, decreasing = T),]
ybr_del_dea <- ybr_del_dea %>% filter(fdr <= 0.05)

ybr_del_dea %>%
  filter(orf_name == 'YBR196C')

ybr_del_dea %>%
  filter(fdr <= 0.05, orf_name %in% gcn4_targets$Target.Systematic.Name) %>%
  nrow()

ybr_del_dea %>% filter(fdr <= 0.05) %>%
  group_by(DE) %>%
  summarise(lfc = median(lfc))

ybr_del_dea %>%
  filter(lfc > 0, orf_name %in% gcn4_targets$Target.Systematic.Name[gcn4_targets$Direction == 'positive'])


write.csv(merge(ybr_del_dea[,c(-5,-6,-7,-8)] %>%
  filter(abs(lfc) > 1, fdr <= 0.05, orf_name %in% gcn4_targets$Target.Systematic.Name),
  all_genes[,c('ORF','GENENAME','DESCRIPTION')], by.x = 'orf_name', by.y = 'ORF'), row.names = F,
  file = 'output/figs/ybr_presentation/ybr_del_dea_gcn4_lfc1.csv')

ybr_del_dea %>%
  filter(orf_name %in% c('YGL195W','YDR283C','YKR026C','YEL009C','YFR009W'))


gcn4.deTable <-  matrix(c(62, 804, 250, 5557),
                  nrow = 2,
                  dimnames = list(DE=c("yes","no"),
                                  GCN4_target=c("yes","no")))
gcn4.deTable
fisher.test(gcn4.deTable)


gcn4lfc1.deTable <-  matrix(c(3, 15, 309, 6346),
                        nrow = 2,
                        dimnames = list(DE=c("yes","no"),
                                        GCN4_target=c("yes","no")))
gcn4lfc1.deTable
fisher.test(gcn4lfc1.deTable)

quantile(coexp.ybr$cor[!is.na(coexp.ybr$orf_name)], c(0.5,0.95,0.99))

coexp.gcn4 <- coexp %>%
  filter(gene1 == 'YEL009C' | gene2 == 'YEL009C')
coexp.gcn4$orf_name[!is.na(coexp.gcn4$gene1) & coexp.gcn4$gene1 == 'YEL009C'] <- coexp.gcn4$gene2[!is.na(coexp.gcn4$gene1) & coexp.gcn4$gene1 == 'YEL009C']
coexp.gcn4$orf_name[!is.na(coexp.gcn4$gene2) & coexp.gcn4$gene2 == 'YEL009C'] <- coexp.gcn4$gene1[!is.na(coexp.gcn4$gene2) & coexp.gcn4$gene2 == 'YEL009C']
coexp.gcn4 <- merge(coexp.gcn4, all_genes, by.x = 'orf_name', by.y = 'ORF', all.x = T)
coexp.gcn4 <- coexp.gcn4[order(coexp.gcn4$cor, decreasing = T),]
row.names(coexp.gcn4) <- NULL

quantile(coexp.gcn4$cor[!is.na(coexp.gcn4$orf_name)], c(0.5,0.95,0.99))

plot.degs <- ybr_del_dea %>%
  ggplot(aes(x = lfc, y = -log10p, col = DE)) +
  geom_point() +
  scale_color_manual(name = '',
                     values = c('Down' = '#303F9F',
                               'Up' = '#FFC107'),
                     na.value = '#757575') +
  coord_cartesian(xlim = c(-3,3)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.margin = margin(0,0,0,0),
        legend.position = 'bottom')
ggsave('output/figs/ybr_presentation/ybr_del_volcano.png', plot.degs,
       height = one.c, width = one.c,
       units = 'mm', dpi = 300,
       bg = 'white')


write.csv(merge(ybr_del_dea %>% filter(abs(lfc) > 1, fdr <= 0.05),
                all_genes, by.x = 'orf_name', by.y = 'ORF'), row.names = F,
          file = 'output/figs/ybr_presentation/ybr_del_dea_lfc1.csv')

ybr_del_dea %>%
  filter(abs(lfc) > 1) %>%
  group_by(DE) %>% count()




#####
pddat.cnts %>%
  group_by(ID, replicate) %>%
  summarise(total_peptide_cnts = sum(peptide_cnts)) %>%
  mutate(ID = paste(ID, replicate, sep = '_')) %>%
  ggplot(aes(x = ID, y = total_peptide_cnts)) +
  stat_summary(geom = 'bar') +
  theme(axis.text = element_text(angle = 90))


pddat.cnts %>%
  filter(peptide_cnts > 0) %>%
  group_by(ID, replicate) %>%
  count() %>%
  mutate(ID = paste(ID, replicate, sep = '_')) %>%
  ggplot(aes(x = ID, y = n)) +
  stat_summary(geom = 'bar')

write.csv(merge(pddat.cnts %>%
                  group_by(ID, replicate) %>%
                  summarise(total_peptide_cnts = sum(peptide_cnts)) %>%
                  mutate(ID = paste(ID, replicate, sep = '_')),
                pddat.cnts %>%
                  filter(peptide_cnts > 0) %>%
                  group_by(ID, replicate) %>%
                  count() %>%
                  mutate(ID = paste(ID, replicate, sep = '_')) ,
                by = c('ID','replicate')),
          file = 'output/figs/ybr_presentation/pull_down_peptide_protein_counts.csv')


##### GSEA for YBR DEL MUTANT
# ybr_del_universe <- 
ybr_del_dea <- dea_results[dea_results$contrast == 'YBR3_BY4741_DEL_vs_YBR3_BY4741_WT',]
ybr_del_dea <- ybr_del_dea[order(ybr_del_dea$lfc, decreasing = T),]
ybr_del_dea <- ybr_del_dea %>% filter(fdr <= 0.05, orf_name %notin% c('YBR115C','YCL018W','YEL021W','YLR303W','YOR202W'))

ybr_del_universe <- ybr_del_dea$lfc
names(ybr_del_universe) <- ybr_del_dea$orf_name  

gsea_ybr_del <- GSEA(ybr_del_universe,
                 TERM2GENE = rbind(orfsets,
                                   data.frame(rbind(overlap_gene.all,
                                                    cbind(set = 'JPN', orf_name = jap_genes),
                                                    cbind(set = 'VMA', orf_name = vma_genes),
                                                    cbind(set = 'GET/SND', orf_name = get_snd_genes)))),
                 eps = 1e-10,
                 nPermSimple = 10000,
                 minGSSize = 3,
                 maxGSSize = 400,
                 verbose = TRUE,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)
plot.gsea.ybr.del <- plot_grid(dotplot(gsea_ybr_del, showCategory=10, split=".sign") + facet_grid(.~.sign),
                           emapplot(enrichplot::pairwise_termsim(gsea_ybr_del), showCategory = 10),
                           cnetplot(gsea_ybr_del, categorySize="pvalue", foldChange=ybr_coexp_universe, showCategory = 5),
                           nrow = 3)
# plot.gsea.ybr.del
# ggsave("output/figs/gsea_results_YBR_DEL.jpg", plot.gsea.ybr.del,
#        height = two.c*3, width = two.c*2, units = 'mm',
#        bg = 'white',
#        dpi = 300)


gsea_ybr_del.res <- data.frame(gsea_ybr_del)
head(gsea_ybr_del.res)
gsea_ybr_del.res <- merge(gsea_ybr_del.res, godetails, by = 'ID', all.x = T)
gsea_ybr_del.res <- gsea_ybr_del.res[order(gsea_ybr_del.res$NES, decreasing = T),]
row.names(gsea_ybr_del.res) <- NULL
head(gsea_ybr_del.res)
write.csv(gsea_ybr_del.res, file = 'output/figs/ybr_presentation/ybr_del_transcriptome_gsea.csv', row.names = F)

##### GO Term Enrichment for COIP Data
temp.deg <- overlap_gene.all %>% data.frame() %>% filter(set == 'COIP_TWO_ALL')
temp.deg <- temp.deg$orf_name
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
data.frame(temp.goe)
coip_two_ha_go <- data.frame(temp.goe)
coip_two_ha_go$GeneRatio <- as.numeric(str_split(coip_two_ha_go$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(coip_two_ha_go$GeneRatio,'/',simplify = T)[,2])
coip_two_ha_go$BgRatio <- as.numeric(str_split(coip_two_ha_go$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(coip_two_ha_go$BgRatio,'/',simplify = T)[,2])
write.csv(coip_two_ha_go, file = 'output/figs/ybr_presentation/pull_down_two_all_GO.csv')

plot.ha.emap <- emapplot(enrichplot::pairwise_termsim(temp.goe), showCategory = 30)
ggsave("output/figs/ybr_presentation/pull_down_two_all_emap.jpg", plot.ha.emap,
       height = 400, width = 400, units = 'mm',
       bg = 'white',
       dpi = 300)

temp.kegg <- enrichKEGG(gene         = temp.deg$ENSEMBL,
                        universe     = all_genes$ENSEMBL,
                        organism     = 'sce',
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff  = 0.05)
data.frame(temp.kegg)
plot.ha.emap <- emapplot(enrichplot::pairwise_termsim(temp.kegg), showCategory = 30)
ggsave("output/figs/ybr_presentation/pull_down_two_all_emap.jpg", plot.ha.emap,
       height = 400, width = 400, units = 'mm',
       bg = 'white',
       dpi = 300)


gsea_coexp2 <- gseGO(ybr_coexp_universe,
                     eps = 1e-10,
                     nPermSimple = 10000,
                     minGSSize = 3,
                     maxGSSize = 800,
                     verbose = TRUE,
                     pAdjustMethod = "BH",
                     OrgDb = org.Sc.sgd.db,
                     keyType       = "ENSEMBL",
                     ont           = "ALL",
                     pvalueCutoff = 0.05)
plot.coexp.emap <- emapplot(enrichplot::pairwise_termsim(gsea_coexp2), showCategory = 30)
ggsave("output/figs/ybr_presentation/ybr_gsea_emap.jpg", plot.coexp.emap,
       height = 400, width = 400, units = 'mm',
       bg = 'white',
       dpi = 300)
