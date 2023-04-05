##### YBR196C-A CO-EXPRESSION ANALYSIS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 11/22/2022

##### INITIALIZE
source(file = 'scripts/initialize.R')
orfs.vma <- c('YOR270C','YGR106C','YMR054W','YHR026W','YLR447C','YEL027W','YKL080W','YOR332W','YPR036W','YPL234C','YBR127C','YGR105W','YEL051W','YGR020C','YMR123W','YHR060W')

coexp <- readRDS('/home/sbp29/R/Data/20221122_coexpression.RDS')
coexp <- coexp %>% filter(num_obs > 100)
head(coexp)

all_genes <- unique(coexp$gene1[!is.na(coexp$gene1)])
all_genes <- bitr(all_genes, fromType = "ORF",
                  toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
                  OrgDb = org.Sc.sgd.db)

ybr_coip <- dbGetQuery(conn, 'select * from YBR_COIP_RESULTS')

##### YBR COEXPRESSION MATRIX
coexp.ybr <- coexp %>%
  filter(gene1 == 'YBR196C-A' | gene2 == 'YBR196C-A')
coexp.ybr$orf_name[!is.na(coexp.ybr$gene1) & coexp.ybr$gene1 == 'YBR196C-A'] <- coexp.ybr$gene2[!is.na(coexp.ybr$gene1) & coexp.ybr$gene1 == 'YBR196C-A']
coexp.ybr$orf_name[!is.na(coexp.ybr$gene2) & coexp.ybr$gene2 == 'YBR196C-A'] <- coexp.ybr$gene1[!is.na(coexp.ybr$gene2) & coexp.ybr$gene2 == 'YBR196C-A']
coexp.ybr <- merge(coexp.ybr, all_genes, by.x = 'orf_name', by.y = 'ORF', all.x = T)
coexp.ybr <- coexp.ybr[order(coexp.ybr$cor, decreasing = T),]
row.names(coexp.ybr) <- NULL
coexp.ybr$rank <- as.numeric(row.names(coexp.ybr))

head(coexp.ybr)

quantile(coexp.ybr$cor, c(0.5,0.95,0.99))

coexp.ybr %>%
  filter(orf_name %in% orfs.vma)

coexp.ybr %>%
  filter(GENENAME %in% c('DPP1','LRO1','ARE2'))

coexp.ybr %>%
  filter(GENENAME %in% c('VPS16','VPS33','VPS41','VAM6','PEP3','PEP5'))

##### GO/KEGG ENRICHMENTS
coexp.annotated.ybr <- coexp.ybr %>% filter(!is.na(orf_name))
paste(coexp.annotated.ybr$orf_name[c(301:400)],collapse = ',')


##### GENE SET ENRICHMENT ANALYSIS
### GO
gene_list <- coexp.annotated.ybr$cor[!is.na(coexp.annotated.ybr$ENSEMBL)]
names(gene_list) <- coexp.annotated.ybr$ENSEMBL[!is.na(coexp.annotated.ybr$ENSEMBL)]

gsea_go <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPermSimple = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             verbose = TRUE, 
             OrgDb = org.Sc.sgd.db,
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05)
head(data.frame(gsea_go))

plot.coexp.gsea.goe <- plot_grid(dotplot(gsea_go, showCategory=10, split=".sign") + facet_grid(.~.sign),
                           emapplot(enrichplot::pairwise_termsim(gsea_go), showCategory = 10),
                           cnetplot(gsea_go, categorySize="pvalue", foldChange=gene_list, showCategory = 5),
                           nrow = 3)

ggsave("output/figs/gsea_go_results_ybr_coexp.jpg", plot.coexp.gsea.goe,
       height = two.c*3, width = two.c*2, units = 'mm',
       bg = 'white',
       dpi = 300)

### KEGG
gsea_kegg <- gseKEGG(geneList=gene_list, 
                nPermSimple = 10000,
                minGSSize = 3,
                maxGSSize = 800,
                verbose = TRUE, 
                organism     = 'sce',
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)
head(data.frame(gsea_kegg))

dotplot(kegg_go, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

xx3 <- enrichplot::pairwise_termsim(kegg_go)
emapplot(xx3, showCategory = 5)

cnetplot(kegg_go, categorySize="pvalue", foldChange=gene_list)


plot.coexp.gsea.kegg <- plot_grid(dotplot(gsea_kegg, showCategory=10, split=".sign") + facet_grid(.~.sign),
                                 emapplot(enrichplot::pairwise_termsim(gsea_kegg), showCategory = 10),
                                 cnetplot(gsea_kegg, categorySize="pvalue", foldChange=gene_list, showCategory = 5),
                                 nrow = 3)

ggsave("output/figs/gsea_kegg_results_ybr_coexp.jpg", plot.coexp.gsea.kegg,
       height = two.c*3, width = two.c*2, units = 'mm',
       bg = 'white',
       dpi = 300)

coexp.gsea_go <- data.frame(gsea_go)
coexp.gsea_go$enrich_direction[coexp.gsea_go$enrichmentScore > 0] <- 'activated'
coexp.gsea_go$enrich_direction[coexp.gsea_go$enrichmentScore < 0] <- 'suppressed'

coexp.gsea_kegg <- data.frame(gsea_kegg)
coexp.gsea_kegg$enrich_direction[coexp.gsea_kegg$enrichmentScore > 0] <- 'activated'
coexp.gsea_kegg$enrich_direction[coexp.gsea_kegg$enrichmentScore < 0] <- 'suppressed'


dbWriteTable(conn, 'GSEA_GO_YBR_COEXPRESSION', coexp.gsea_go, overwrite = T)
dbWriteTable(conn, 'GSEA_KEGG_YBR_COEXPRESSION', coexp.gsea_kegg, overwrite = T)

# ##### YBR and VATPases
# vma1_interactors <- read.table(file = 'input/VMA1_physical_interactions.txt', header = T, sep = '\t', skip = 8)
# ybr_coip$gene[ybr_coip$gene %in% unique(vma1_interactors$Interactor.1)]
# 
# vma6_interactors <- read.table(file = 'input/VMA6_physical_interactions.txt', header = T, sep = '\t', skip = 8)
# ybr_coip$gene[ybr_coip$gene %in% unique(vma6_interactors$Interactor.1)]
# 
# ybr_coip[ybr_coip$gene %in% c('HSC82','RVB1'),]


##### COIP and GSEA
### GO
coip_gene_list <- ybr_coip$score
names(coip_gene_list) <- ybr_coip$ORF

gsea_go <- gseGO(geneList=coip_gene_list, 
                ont ="ALL", 
                keyType = "ENSEMBL", 
                nPermSimple = 10000,
                minGSSize = 3,
                maxGSSize = 100,
                verbose = TRUE, 
                OrgDb = org.Sc.sgd.db,
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)

plot.coip.gsea.goe <- plot_grid(dotplot(gsea_go, showCategory=10, split=".sign") + facet_grid(.~.sign),
                                 emapplot(enrichplot::pairwise_termsim(gsea_go), showCategory = 10),
                                 cnetplot(gsea_go, categorySize="pvalue", foldChange=gene_list, showCategory = 5),
                                 nrow = 3)
ggsave("output/figs/gsea_go_results_ybr_coip.jpg", plot.coip.gsea.goe,
       height = two.c*3, width = two.c*2, units = 'mm',
       bg = 'white',
       dpi = 300)

### KEGG
gsea_kegg <- gseKEGG(geneList=coip_gene_list, 
                   nPermSimple = 10000,
                   minGSSize = 3,
                   maxGSSize = 100,
                   verbose = TRUE, 
                   organism     = 'sce',
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)
head(data.frame(gsea_kegg))


coip.gsea_go <- data.frame(gsea_go)
coip.gsea_go$enrich_direction[coip.gsea_go$enrichmentScore > 0] <- 'activated'
coip.gsea_go$enrich_direction[coip.gsea_go$enrichmentScore < 0] <- 'suppressed'

dbWriteTable(conn, 'GSEA_GO_YBR_COIP', coip.gsea_go, overwrite = T)



##### EMPIRICAL SET ENRICHMENT
gene_list <- coexp.annotated.ybr$cor[!is.na(coexp.annotated.ybr$ENSEMBL)]
names(gene_list) <- coexp.annotated.ybr$ENSEMBL[!is.na(coexp.annotated.ybr$ENSEMBL)]

# hello <- GSEA(gene_list,
#               TERM2GENE = data.frame(rbind(cbind(term = 'COIP', name = ybr_coip$ORF),
#                                            cbind(term = 'VMA', name = orfs.vma))),
#               eps = 1e-10,
#               nPermSimple = 10000,
#               # minGSSize = 3,
#               # maxGSSize = 200,
#               verbose = TRUE, 
#               pAdjustMethod = "BH",
#               pvalueCutoff = 0.05)
# data.frame(hello)

go2orf <- as.list(org.Sc.sgdGO2ORF)
orfsets <- NULL
for (i in 1:length(go2orf)) {
  temp <- data.frame(set = names(go2orf[i]),
                     genes = melt(go2orf[[i]]))
  orfsets <- rbind(orfsets, temp)
}
colnames(orfsets) <- c('set','orf_name')
head(orfsets)

GO <- as.list(GOTERM)
godetails <- NULL
for (i in 1:length(GO)) {
  godetails <- rbind(godetails,
                     data.frame(ID = GOID(GO[[i]]),
                           ONT = Ontology(GO[[i]]),
                           Term = Term(GO[[i]]),
                           Definition = Definition(GO[[i]])))
}
head(godetails)


gsea_ybr <- GSEA(gene_list,
              TERM2GENE = rbind(orfsets,
                                data.frame(rbind(cbind(set = 'COIP', orf_name = ybr_coip$ORF),
                                                 cbind(set = 'VMA', orf_name = orfs.vma),
                                                 cbind(set = 'GET/SND', orf_name = c('YGL020C','YER083C','YDL100C',
                                                                                     'YLR065C','YBR106W'))))),
              eps = 1e-10,
              nPermSimple = 10000,
              minGSSize = 3,
              maxGSSize = 800,
              verbose = TRUE,
              pAdjustMethod = "BH",
              pvalueCutoff = 0.05)
plot.gsea.ybr <- plot_grid(dotplot(gsea_ybr, showCategory=10, split=".sign") + facet_grid(.~.sign),
                            emapplot(enrichplot::pairwise_termsim(gsea_ybr), showCategory = 10),
                            cnetplot(gsea_ybr, categorySize="pvalue", foldChange=gene_list, showCategory = 5),
                            nrow = 3)
ggsave("output/figs/gsea_kegg_results_YBR_COEXP.jpg", plot.gsea.ybr,
       height = two.c*3, width = two.c*2, units = 'mm',
       bg = 'white',
       dpi = 300)

gsea_ybr <- data.frame(gsea_ybr)
rownames(gsea_ybr) <- NULL
gsea_ybr <- merge(gsea_ybr, godetails, by = 'ID', all.x = T)
gsea_ybr <- gsea_ybr[order(gsea_ybr$NES, gsea_ybr$p.adjust, gsea_ybr$qvalues, decreasing = T),]
rownames(gsea_ybr) <- NULL

dbWriteTable(conn, 'GSEA_YBR_COEXP', gsea_ybr, overwrite = T)


#####
gsea_ybr[str_detect(gsea_ybr$core_enrichment, 'YMR013C'),]
coexp.ybr %>%
  filter(orf_name == 'YMR013C')

str_replace_all(gsea_ybr$core_enrichment[gsea_ybr$Term == 'cation transmembrane transport' & !is.na(gsea_ybr$Term)],
            '/', ', ')

str_split(gsea_ybr$core_enrichment[gsea_ybr$Term == 'endoplasmic reticulum membrane' & !is.na(gsea_ybr$Term)], '/')

gsea_ybr %>%
  filter(ONT == 'BP')


#####
coexp.ybr %>%
  filter(orf_name %in% str_split(gsea_ybr$core_enrichment[gsea_ybr$Term == 'endoplasmic reticulum membrane' & !is.na(gsea_ybr$Term)], '/', simplify = T)) %>%
  ggplot() +
  geom_vline(aes(xintercept = rank)) +
  coord_cartesian(xlim = c(1,max(coexp.ybr$rank)))

str_split(gsea_ybr$core_enrichment[gsea_ybr$ID == 'COIP'], '/', simplify = T)

coexp.ybr %>%
  filter(orf_name %in% str_split(gsea_ybr$core_enrichment[gsea_ybr$ID == 'COIP'], '/', simplify = T)) %>%
  ggplot() +
  geom_vline(aes(xintercept = rank)) +
  coord_cartesian(xlim = c(1,max(coexp.ybr$rank)))

coexp.ybr %>%
  filter(orf_name %in% ybr_coip$ORF) %>%
  ggplot() +
  geom_vline(aes(xintercept = rank)) +
  coord_cartesian(xlim = c(1,max(coexp.ybr$rank)))

