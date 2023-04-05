##### YBR196C-A ALL SHORT READ TRANSCRIPTOME
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 11/15/2022

source(file = 'scripts/initialize.R')

##### GATHER DATA
dea_results <- dbGetQuery(conn, "select * from RNASEQ_YBR_DEA a where 
        ((a.batch = 'YBR1' and a.comparison in ('ATG_vs_WT', 'DEL_vs_WT')) or
        (a.batch = 'YBR2' and a.comparison in ('DEL_vs_WT')) or
        (a.batch = 'YBR3' and a.comparison in ('ATG_vs_WT', 'DEL_vs_WT')) or
        (a.batch = 'YBR1/3' and a.comparison in ('ATG_vs_WT')))")

##### LFC Ranks
for (i in unique(dea_results$contrast)) {
  dea_results$rank[dea_results$contrast == i] <- order(dea_results$lfc[dea_results$contrast == i], decreasing = T)
}

##### LFC value and rank correlation
head(dea_results)
contrast_cmbn <- t(combn(unique(dea_results$contrast), 2))

dea_correlation <- NULL
for (i in seq(1,nrow(contrast_cmbn))) {
  temp <- merge(dea_results[dea_results$contrast == contrast_cmbn[i,1],],
                dea_results[dea_results$contrast == contrast_cmbn[i,2],],
                by = 'orf_name', suffixes = c('_set1','_set2'))
  dea_correlation <- rbind(dea_correlation, temp)
}
head(dea_correlation)

dea_correlation %>%
  ggplot(aes(x = rank_set1, y = rank_set2)) +
  geom_point(size = 0.5) +
  geom_density_2d() +
  facet_wrap(contrast_set1 ~ contrast_set2) 


plot.lfc.corr <- dea_correlation %>%
  ggplot(aes(x = lfc_set1, y = lfc_set2)) +
  geom_point(size = 0.5) +
  geom_point(data = dea_correlation %>%
               filter(!is.na(DE_set1), !is.na(DE_set2)),
             aes(x = lfc_set1, y = lfc_set2), size = 0.5, col = 'red') +
  facet_wrap(batch_set1*background_set1*comparison_set1 ~ batch_set2*background_set2*comparison_set2) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        strip.text = element_text(size = txt-2,
                                  margin = margin(1,0,0.1,0, "mm")))
ggsave(sprintf("output/figs/rnaseq_ybr_sr_de_lfc_corr.jpg"), plot.lfc.corr,
       height = two.c, width = two.c*3/2, units = 'mm',
       bg = 'white',
       dpi = 300)

##### GO KEGG ENRICHMENT IN DEGS
all_genes <- unique(dea_results$orf_name)
all_genes <- bitr(all_genes, fromType = "ORF",
                  toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
                  OrgDb = org.Sc.sgd.db)

head(dea_results)
goe <- NULL
kegg <- NULL
for (c in unique(dea_results$contrast)) {
  for (de in unique(dea_results$DE[dea_results$contrast == c & !is.na(dea_results$DE)])) {
    batch <- unique(dea_results$batch[dea_results$contrast == c & dea_results$DE == de & !is.na(dea_results$DE)])
    bg <- unique(dea_results$background[dea_results$contrast == c & dea_results$DE == de & !is.na(dea_results$DE)])
    com <- unique(dea_results$comparison[dea_results$contrast == c & dea_results$DE == de & !is.na(dea_results$DE)])
    
    temp.deg <- dea_results$orf_name[dea_results$contrast == c & dea_results$DE == de & !is.na(dea_results$DE)]
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
      cat(sprintf('There are no GO term enrichment for %s DEGs in %s.\n\n',de,c))
    } else{
      cat(sprintf('Top 5 GO term enrichment of %s DEGs in %s for:\n',de,c))
      cat(sprintf('\t1. Cellular component - %s.\n',
                  paste(head(unique(temp.goe$Description[temp.goe$ONTOLOGY == 'CC']),5),collapse = ', ')))
      cat(sprintf('\t2. Biological process - %s.\n',
                  paste(head(unique(temp.goe$Description[temp.goe$ONTOLOGY == 'BP']),5),collapse = ', ')))
      cat(sprintf('\t3. Molecular function - %s.\n\n',
                  paste(head(unique(temp.goe$Description[temp.goe$ONTOLOGY == 'MF']),5),collapse = ', ')))
      goe <- rbind(goe, data.frame(temp.goe, contrast = c, batch = batch, background = bg, comparison = com, DE = de))
    }
    
    temp.kegg <- enrichKEGG(gene         = temp.deg$ENSEMBL,
                            universe     = all_genes$ENSEMBL,
                            organism     = 'sce',
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff  = 0.05)
    if (dim(temp.kegg)[1] == 0) {
      cat(sprintf('There are no KEGG pathway enriched for %s DEGs in %s.\n\n',de,c))
    } else{
      cat(sprintf('Top 10 KEGG pathway enrichment for %s DEGs in %s are:\n%s.\n\n',
                  de,c,
                  paste(head(unique(temp.kegg$Description),10),collapse = ', ')))
      kegg <- rbind(kegg, data.frame(temp.kegg, contrast = c, batch = batch, background = bg, comparison = com, DE = de))
    }
  }
}
goe$GeneRatio <- as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,2])
goe$BgRatio <- as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,2])
goe$GO <- paste0(goe$ONTOLOGY, '_', goe$Description)
goe <- goe[order(goe$contrast, goe$batch, goe$background, goe$comparison, goe$DE, goe$GeneRatio, goe$qvalue, decreasing = T),]

kegg$GeneRatio <- as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,2])
kegg$BgRatio <- as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,2])
kegg <- kegg[order(kegg$contrast, kegg$batch, kegg$background, kegg$comparison, kegg$DE, kegg$GeneRatio, kegg$qvalue, decreasing = T),]

dbWriteTable(conn, 'RNASEQ_YBR_GOE', goe, overwrite = T)
dbWriteTable(conn, 'RNASEQ_YBR_KEGG', kegg, overwrite = T)


head(goe)

goe %>%
  group_by(DE,GO) %>%
  count() %>%
  filter(n > 2)

kegg %>%
  group_by(DE, Description) %>%
  count() %>%
  filter(n > 2)

paste(str_replace_all(kegg$geneID[kegg$Description == 'Carbon metabolism'], '/',','), collapse = ',')
hello <- melt(str_split(kegg$geneID[kegg$Description == 'Carbon metabolism'], '/', simplify = T))
hello %>%
  group_by(value) %>%
  count() %>%
  filter(n > 2)

dea_results %>%
  filter(orf_name == 'YBR196C')

dea_common <- dea_results %>%
  filter(!is.na(DE)) %>%
  group_by(DE, orf_name) %>%
  count() %>% filter(n > 3) %>% 
  data.frame()
temp <- bitr(dea_common$orf_name, fromType = "ORF",
             toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)

# ##### OVERLAP USING UPSET PLOTS
# common_up = list("YBR1 BY4742 DEL vs WT" = dea_results$orf_name[!is.na(dea_results$DE) &dea_results$DE == 'Up' &
#                                                                   dea_results$contrast == "YBR1_BY4742_DEL_vs_YBR1_BY4742_WT"],
#                  "YBR1 BY4742 ATG vs WT" = dea_results$orf_name[!is.na(dea_results$DE) &dea_results$DE == 'Up' &
#                                                                   dea_results$contrast == "YBR1_BY4742_ATG_vs_YBR1_BY4742_WT"],
#                  "YBR3 BY4742 ATG vs WT" = dea_results$orf_name[!is.na(dea_results$DE) &dea_results$DE == 'Up' &
#                                                                   dea_results$contrast == "YBR3_BY4742_ATG_vs_YBR3_BY4742_WT"],
#                  "YBR2 BY4742 DEL vs WT" = dea_results$orf_name[!is.na(dea_results$DE) &dea_results$DE == 'Up' &
#                                                                   dea_results$contrast == "YBR2_BY4742_DEL_vs_YBR2_BY4742_WT"],
#                  "YBR3 BY4741 DEL vs WT" = dea_results$orf_name[!is.na(dea_results$DE) &dea_results$DE == 'Up' &
#                                                                   dea_results$contrast == "YBR3_BY4741_DEL_vs_YBR3_BY4741_WT"],
#                  "YBR1/3 BY4742 ATG vs WT" = dea_results$orf_name[!is.na(dea_results$DE) &dea_results$DE == 'Up' &
#                                                                   dea_results$contrast == "YBR1_YBR3_BY4742_ATG_vs_YBR1_YBR3_BY4742_WT"])
# 
# upset(fromList(common_up), order.by = "freq",nsets = 6)
# 
# 
# common_down = list("YBR1 BY4742 DEL vs WT" = dea_results$orf_name[!is.na(dea_results$DE) &dea_results$DE == 'Down' &
#                                                                   dea_results$contrast == "YBR1_BY4742_DEL_vs_YBR1_BY4742_WT"],
#                  "YBR1 BY4742 ATG vs WT" = dea_results$orf_name[!is.na(dea_results$DE) &dea_results$DE == 'Down' &
#                                                                   dea_results$contrast == "YBR1_BY4742_ATG_vs_YBR1_BY4742_WT"],
#                  "YBR3 BY4742 ATG vs WT" = dea_results$orf_name[!is.na(dea_results$DE) & dea_results$DE == 'Down' &
#                                                                   dea_results$contrast == "YBR3_BY4742_ATG_vs_YBR3_BY4742_WT"],
#                  "YBR2 BY4742 DEL vs WT" = dea_results$orf_name[!is.na(dea_results$DE) &dea_results$DE == 'Down' &
#                                                                   dea_results$contrast == "YBR2_BY4742_DEL_vs_YBR2_BY4742_WT"],
#                  "YBR3 BY4741 DEL vs WT" = dea_results$orf_name[!is.na(dea_results$DE) &dea_results$DE == 'Down' &
#                                                                   dea_results$contrast == "YBR3_BY4741_DEL_vs_YBR3_BY4741_WT"],
#                  "YBR1/3 BY4742 ATG vs WT" = dea_results$orf_name[!is.na(dea_results$DE) &dea_results$DE == 'Down' &
#                                                                     dea_results$contrast == "YBR1_YBR3_BY4742_ATG_vs_YBR1_YBR3_BY4742_WT"])
# 
# upset(fromList(common_down), order.by = "freq",nsets = 6)
# 
# dea_results %>%
#   filter(!is.na(DE)) %>%
#   group_by(DE, orf_name) %>%
#   count() %>%
#   filter(n > 4)
# 
# dea_results %>%
#   filter(orf_name == 'YBR196C')


##### GSEA Analysis
gsea <- NULL
for (c in unique(dea_results$contrast)) {
  gene_list <- dea_results$lfc[dea_results$contrast == c & dea_results$fdr <= 0.05] # 
  names(gene_list) <- dea_results$orf_name[dea_results$contrast == c & dea_results$fdr <= 0.05]
  
  gse <- gseGO(geneList=gene_list, 
               ont ="ALL", 
               keyType = "ENSEMBL", 
               nPermSimple = 10000,
               minGSSize = 3,
               maxGSSize = 400,
               verbose = TRUE, 
               OrgDb = org.Sc.sgd.db,
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05)
  head(data.frame(gse))
  
  dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
  }



