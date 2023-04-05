##### YBR196C-A ALL SHORT READ TRANSCRIPTOME
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 12/01/2022

source(file = 'scripts/initialize.R')

##### GATHER DATA
dea_results <- dbGetQuery(conn, "select * from RNASEQ_YBR_DEA a where 
        ((a.batch = 'YBR1' and a.comparison in ('ATG_vs_WT', 'DEL_vs_WT')) or
        (a.batch = 'YBR2' and a.comparison in ('ATG_vs_WT', 'DEL_vs_WT')) or
        (a.batch = 'YBR3' and a.comparison in ('ATG_vs_WT', 'DEL_vs_WT')))")

dea_results %>%
  filter(!is.na(DE), abs(lfc) > 1) %>%
  group_by(contrast, DE, orf_name) %>%
  count() %>% data.frame()


##### GSEA Analysis
## GSEA GO
all_genes <- unique(dea_results$orf_name)
all_genes <- bitr(all_genes, fromType = "ORF",
                  toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
                  OrgDb = org.Sc.sgd.db)

dea_results <- merge(dea_results, all_genes, by.x = 'orf_name', by.y = 'ORF', all.x = T)
dea_results <- dea_results[order(dea_results$contrast, dea_results$lfc, decreasing = T),]
head(dea_results)

gsea_go <- NULL
gsea_kegg <- NULL
for (c in unique(dea_results$contrast)) {
  batch <- unique(dea_results$batch[dea_results$contrast == c])
  bg <- unique(dea_results$background[dea_results$contrast == c])
  com <- unique(dea_results$comparison[dea_results$contrast == c])
  
  for (de in unique(dea_results$DE[dea_results$contrast == c & !is.na(dea_results$DE)])) {
    temp <- dea_results[!is.na(dea_results$ENSEMBL) & dea_results$contrast == c & dea_results$DE == de & !is.na(dea_results$DE),]
    temp <- temp[order(abs(temp$lfc), decreasing = T),]
    gene_list <- abs(temp$lfc[!is.na(temp$ENSEMBL)])
    names(gene_list) <- temp$orf_name[!is.na(temp$ENSEMBL)]

    temp.goe <- gseGO(geneList=gene_list,
                      ont ="ALL",
                      keyType = "ENSEMBL",
                      nPermSimple = 10000,
                      minGSSize = 3,
                      maxGSSize = 400,
                      verbose = TRUE,
                      OrgDb = org.Sc.sgd.db,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05)

    if (dim(data.frame(temp.goe))[1] > 0) {
      # plot.temp.goe <- plot_grid(dotplot(temp.goe, showCategory=10, split=".sign") + facet_grid(.~.sign),
      #                            emapplot(enrichplot::pairwise_termsim(temp.goe), showCategory = 10),
      #                            cnetplot(temp.goe, categorySize="pvalue", foldChange=gene_list, showCategory = 5),
      #                            nrow = 3)
      # ggsave(sprintf("output/figs/gsea_go_results_%s_%s.jpg",de,c), plot.temp.goe,
      #        height = two.c*3, width = two.c*2, units = 'mm',
      #        bg = 'white',
      #        dpi = 300)

      gsea_go <- rbind(gsea_go, cbind(batch = batch, background = bg, comparison = com, contrast = c, DE = de, data.frame(temp.goe)))
    }


    temp.kegg <- gseKEGG(geneList=gene_list,
                         nPermSimple = 10000,
                         minGSSize = 3,
                         maxGSSize = 400,
                         verbose = TRUE,
                         organism     = 'sce',
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05)
    head(data.frame(temp.kegg))

    if (dim(data.frame(temp.kegg))[1] > 0) {
      # plot.temp.kegg <- plot_grid(dotplot(temp.kegg, showCategory=10, split=".sign") + facet_grid(.~.sign),
      #                             emapplot(enrichplot::pairwise_termsim(temp.kegg), showCategory = 10),
      #                             cnetplot(temp.kegg, categorySize="pvalue", foldChange=gene_list, showCategory = 5),
      #                             nrow = 3)
      # ggsave(sprintf("output/figs/gsea_kegg_results_%s_%s.jpg",de,c), plot.temp.kegg,
      #        height = two.c*3, width = two.c*2, units = 'mm',
      #        bg = 'white',
      #        dpi = 300)

      gsea_kegg <- rbind(gsea_kegg, cbind(batch = batch, background = bg, comparison = com, contrast = c, DE = de, data.frame(temp.kegg)))
    }
  }
}

gsea_go$enrichmentScore[gsea_go$DE == 'Down'] <- gsea_go$enrichmentScore[gsea_go$DE == 'Down'] * -1
gsea_go$enrich_direction[gsea_go$enrichmentScore > 0] <- 'activated'
gsea_go$enrich_direction[gsea_go$enrichmentScore < 0] <- 'suppressed'

gsea_kegg$enrichmentScore[gsea_kegg$DE == 'Down'] <- gsea_kegg$enrichmentScore[gsea_kegg$DE == 'Down'] * -1
gsea_kegg$enrich_direction[gsea_kegg$enrichmentScore > 0] <- 'activated'
gsea_kegg$enrich_direction[gsea_kegg$enrichmentScore < 0] <- 'suppressed'


dbWriteTable(conn, 'GSEA_GO_YBR_RNASEQ', gsea_go, overwrite = T)
dbWriteTable(conn, 'GSEA_KEGG_YBR_RNASEQ', gsea_kegg, overwrite = T)


