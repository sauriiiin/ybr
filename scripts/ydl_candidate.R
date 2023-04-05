coexp <- readRDS('/home/sbp29/R/Data/20221122_coexpression.RDS')
coexp.ydl <- coexp %>%
  filter(row == 'chr4_94133' | column == 'chr4_94133')
head(coexp.ydl)
coexp.ydl <- coexp.ydl[!is.na(coexp.ydl$gene2),]
coexp.ydl$orf_name <- coexp.ydl$gene2

coexp.ydl <- merge(coexp.ydl, all_genes, by.x = 'orf_name', by.y = 'ORF', all.x = T)
coexp.ydl <- coexp.ydl[order(coexp.ydl$cor, decreasing = T),]
row.names(coexp.ydl) <- NULL

head(coexp.ydl)


ydl_coexp_universe <- coexp.ydl$cor[!is.na(coexp.ydl$ENSEMBL)]
names(ydl_coexp_universe) <- coexp.ydl$ENSEMBL[!is.na(coexp.ydl$ENSEMBL)]


gsea_ydl <- gseGO(ydl_coexp_universe,
                    eps = 0,
                    nPermSimple = 10000,
                    minGSSize = 3,
                    maxGSSize = 800,
                    verbose = TRUE,
                    pAdjustMethod = "BH",
                    # organism     = 'sce',
                    OrgDb = org.Sc.sgd.db,
                    keyType       = "ENSEMBL",
                    ont           = "ALL",
                    pvalueCutoff = 0.05)


gsea_ydl.res <- data.frame(gsea_ydl)
plot.coexp.emap <- emapplot(enrichplot::pairwise_termsim(gsea_ydl), showCategory = 30)
ggsave("output/figs/ydl_gsea_emap.jpg", plot.coexp.emap,
       height = 400, width = 400, units = 'mm',
       bg = 'white',
       dpi = 300)
?emapplot
data.frame(gsea_ydl)

dotplot(gsea_ydl, showCategory=10, split=".sign") + facet_grid(.~.sign)
plot.coexp.tree <- treeplot(enrichplot::pairwise_termsim(gsea_ydl), hclust_method = "average")
ggsave("output/figs/ydl_gsea_tree.jpg", plot.coexp.tree,
       height = 200, width = 400, units = 'mm',
       bg = 'white',
       dpi = 300)


gsea_ydl.res <- gsea_ydl.res[order(gsea_ydl.res$NES, decreasing = T),]

write.csv(coexp.ydl, file = 'output/ydl_coexpression_network.csv')
write.csv(gsea_ydl.res, file = 'output/ydl_coexpression_gsea_kegg.csv')
