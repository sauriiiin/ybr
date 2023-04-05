
library(RMariaDB)
library(cowplot)
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")

conn <- initialize.sql("saurin_test")

all_genes <- dbGetQuery(conn, 'select * from RNASEQ_YBRATG_ALL_GENES')
vma_genes <- c('YOR270C','YGR106C','YMR054W','YHR026W','YLR447C','YEL027W','YKL080W','YOR332W','YPR036W','YPL234C','YBR127C','YGR105W','YEL051W','YGR020C','YMR123W','YHR060W')


source('/home/sbp29/R/Projects/ybr/scripts/pull_down_results/analysis.R')
load(file = '/home/sbp29/R/Projects/ybr/scripts/pull_down_results/pddata.RData')
unique(pddat$ID)

combine_replicates <- 'Y'

test_group <- c("attempt_2_WT_YBR_mNG", "attempt_2_doa10D_YBR_mNG", "attempt_1_doa10D_YBR_mNG", "attempt_2_WT_YBR_HA")
control_group <- c("attempt_2_WT_mNG", "attempt_2_doa10D_mNG", "attempt_1_WT_mNG", "attempt_2_WT_HA" )

resres <- ms_overlap(combine_replicates, test_group, control_group)

upset_plot <- resres[[1]]
result_summary <- resres[[2]]
result_proteins <- resres[[3]]
results_go <- resres[[4]]
results_kegg <- resres[[5]]


overlap_genes <- result_proteins %>%
  filter(overlap == 1) %>%
  group_by(Gene) %>%
  count() %>% data.frame()


overlap_gene.descriptions <- all_genes[(all_genes$ORF %in% overlap_genes$Gene) |
            (all_genes$GENENAME %in% overlap_genes$Gene),]


sum((all_genes$ORF %in% overlap_genes$Gene) |
      (all_genes$GENENAME %in% overlap_genes$Gene))


##### COEXP DATA
coexp <- readRDS('/home/sbp29/R/Data/20221122_coexpression.RDS')
coexp <- coexp %>% filter(num_obs > 100)
head(coexp)

coexp_universe <- unique(coexp$gene1[!is.na(coexp$gene1)])
coexp_universe <- bitr(coexp_universe, fromType = "ORF",
                  toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
                  OrgDb = org.Sc.sgd.db)

coexp.ybr <- coexp %>%
  filter(gene1 == 'YBR196C-A' | gene2 == 'YBR196C-A')
coexp.ybr$orf_name[!is.na(coexp.ybr$gene1) & coexp.ybr$gene1 == 'YBR196C-A'] <- coexp.ybr$gene2[!is.na(coexp.ybr$gene1) & coexp.ybr$gene1 == 'YBR196C-A']
coexp.ybr$orf_name[!is.na(coexp.ybr$gene2) & coexp.ybr$gene2 == 'YBR196C-A'] <- coexp.ybr$gene1[!is.na(coexp.ybr$gene2) & coexp.ybr$gene2 == 'YBR196C-A']
coexp.ybr <- merge(coexp.ybr, all_genes, by.x = 'orf_name', by.y = 'ORF', all.x = T)
coexp.ybr <- coexp.ybr[order(coexp.ybr$cor, decreasing = T),]
row.names(coexp.ybr) <- NULL

quantile(coexp.ybr$cor, c(0.5,0.95,0.99))

##### COEXPRESSION RESULTS FOR THE OVERLAP GENES
coexp.ybr %>%
  filter(orf_name %in% overlap_gene.descriptions$ORF)

##### GSEA WITH COEXP OF THE OVERLAP RESULTS
# check ybr_coep.R script

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

coexp.annotated.ybr <- coexp.ybr %>% filter(!is.na(orf_name))
ybr_coexp_universe <- coexp.annotated.ybr$cor[!is.na(coexp.annotated.ybr$ENSEMBL)]
names(ybr_coexp_universe) <- coexp.annotated.ybr$ENSEMBL[!is.na(coexp.annotated.ybr$ENSEMBL)]

gsea_ybr <- GSEA(ybr_coexp_universe,
                 TERM2GENE = rbind(orfsets,
                                   data.frame(rbind(
                                     cbind(set = 'COIP', orf_name = overlap_gene.descriptions$ORF),
                                     cbind(set = 'VMA', orf_name = vma_genes),
                                     cbind(set = 'GET/SND', orf_name = c('YGL020C','YER083C','YDL100C',
                                                                         'YLR065C','YBR106W'))))),
                 eps = 1e-10,
                 nPermSimple = 10000,
                 minGSSize = 3,
                 maxGSSize = 800,
                 verbose = TRUE,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)
# plot.gsea.ybr <- plot_grid(dotplot(gsea_ybr, showCategory=10, split=".sign") + facet_grid(.~.sign),
#                            emapplot(enrichplot::pairwise_termsim(gsea_ybr), showCategory = 10),
#                            cnetplot(gsea_ybr, categorySize="pvalue", foldChange=ybr_coexp_universe, showCategory = 5),
#                            nrow = 3)
# plot.gsea.ybr 
# ggsave("output/figs/gsea_kegg_results_YBR_COEXP.jpg", plot.gsea.ybr,
#        height = two.c*3, width = two.c*2, units = 'mm',
#        bg = 'white',
#        dpi = 300)

gsea_ybr <- data.frame(gsea_ybr)
rownames(gsea_ybr) <- NULL
gsea_ybr <- merge(gsea_ybr, godetails, by = 'ID', all.x = T)
gsea_ybr <- gsea_ybr[order(gsea_ybr$NES, gsea_ybr$p.adjust, gsea_ybr$qvalues, decreasing = T),]
rownames(gsea_ybr) <- NULL

