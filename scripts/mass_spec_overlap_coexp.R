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
`%notin%` <- Negate(`%in%`)
# source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
# source('/home/sbp29/R/Projects/ybr/scripts/pull_down_results/analysis.R')
# load(file = '/home/sbp29/R/Projects/ybr/scripts/pull_down_results/pddata.RData')
# coexp <- readRDS('/home/sbp29/R/Data/20221122_coexpression.RDS')
# conn <- initialize.sql("saurin_test")
# 
# all_genes <- dbGetQuery(conn, 'select * from RNASEQ_YBRATG_ALL_GENES')
# vma_genes <- c('YOR270C','YGR106C','YMR054W','YHR026W','YLR447C','YEL027W','YKL080W','YOR332W',
#                'YPR036W','YPL234C','YBR127C','YGR105W','YEL051W','YGR020C','YMR123W','YHR060W')
# jap_genes <- c('YBR109C','YJR106W','YBR187W')
# get_snd_genes <- c('YGL020C','YER083C','YDL100C','YLR065C','YBR106W')
# 
# ##### PROCESS COEXP DATA
# coexp <- coexp %>% filter(num_obs > 100)
# head(coexp)
# 
# coexp_universe <- unique(coexp$gene1[!is.na(coexp$gene1)])
# coexp_universe <- bitr(coexp_universe, fromType = "ORF",
#                        toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
#                        OrgDb = org.Sc.sgd.db)
# 
# coexp.ybr <- coexp %>%
#   filter(gene1 == 'YBR196C-A' | gene2 == 'YBR196C-A')
# coexp.ybr$orf_name[!is.na(coexp.ybr$gene1) & coexp.ybr$gene1 == 'YBR196C-A'] <- coexp.ybr$gene2[!is.na(coexp.ybr$gene1) & coexp.ybr$gene1 == 'YBR196C-A']
# coexp.ybr$orf_name[!is.na(coexp.ybr$gene2) & coexp.ybr$gene2 == 'YBR196C-A'] <- coexp.ybr$gene1[!is.na(coexp.ybr$gene2) & coexp.ybr$gene2 == 'YBR196C-A']
# coexp.ybr <- merge(coexp.ybr, all_genes, by.x = 'orf_name', by.y = 'ORF', all.x = T)
# coexp.ybr <- coexp.ybr[order(coexp.ybr$cor, decreasing = T),]
# row.names(coexp.ybr) <- NULL
# 
# # quantile(coexp.ybr$cor, c(0.5,0.95,0.99))
# 
# ##### COMPILE GO DATA
# go2orf <- as.list(org.Sc.sgdGO2ORF)
# orfsets <- NULL
# for (i in 1:length(go2orf)) {
#   temp <- data.frame(set = names(go2orf[i]),
#                      genes = melt(go2orf[[i]]))
#   orfsets <- rbind(orfsets, temp)
# }
# colnames(orfsets) <- c('set','orf_name')
# head(orfsets)
# 
# GO <- as.list(GOTERM)
# godetails <- NULL
# for (i in 1:length(GO)) {
#   godetails <- rbind(godetails,
#                      data.frame(ID = GOID(GO[[i]]),
#                                 ONT = Ontology(GO[[i]]),
#                                 Term = Term(GO[[i]]),
#                                 Definition = Definition(GO[[i]])))
# }
# head(godetails)

##### LOAD PRESAVED DATA
load('~/R/Projects/ybr/data/230119_pulldown_rnaseq_coexp_init.RData')
conn <- initialize.sql("saurin_test")

##### PROCESS EXPRESSION DATA
dea_results <- dbGetQuery(conn, "select * from RNASEQ_YBR_DEA a where 
        ((a.batch = 'YBR1' and a.comparison in ('ATG_vs_WT', 'DEL_vs_WT')) or
        (a.batch = 'YBR2' and a.comparison in ('ATG_vs_WT', 'DEL_vs_WT')) or
        (a.batch = 'YBR3' and a.comparison in ('ATG_vs_WT', 'DEL_vs_WT')))")

dea_results.lfc1 <- dea_results %>%
  filter(!is.na(DE), abs(lfc) > 1) %>%
  group_by(contrast, comparison, DE, orf_name) %>%
  count() %>% data.frame()
unique(dea_results.lfc1$contrast)

dea_results.lfc1.orf_sets <- rbind(
  cbind(set = "YBR1_BY4742_ATG_vs_YBR1_BY4742_WT_DOWN",
        orf_name = dea_results.lfc1$orf_name[dea_results.lfc1$contrast == "YBR1_BY4742_ATG_vs_YBR1_BY4742_WT" &
                                               dea_results.lfc1$DE == "Down"]),
  cbind(set = "YBR1_BY4742_ATG_vs_YBR1_BY4742_WT_UP",
        orf_name = dea_results.lfc1$orf_name[dea_results.lfc1$contrast == "YBR1_BY4742_ATG_vs_YBR1_BY4742_WT" &
                                               dea_results.lfc1$DE == "Up"]),
  cbind(set = "YBR1_BY4742_DEL_vs_YBR1_BY4742_WT_DOWN",
        orf_name = dea_results.lfc1$orf_name[dea_results.lfc1$contrast == "YBR1_BY4742_DEL_vs_YBR1_BY4742_WT" &
                                               dea_results.lfc1$DE == "Down"]),
  cbind(set = "YBR1_BY4742_DEL_vs_YBR1_BY4742_WT_UP",
        orf_name = dea_results.lfc1$orf_name[dea_results.lfc1$contrast == "YBR1_BY4742_DEL_vs_YBR1_BY4742_WT" &
                                               dea_results.lfc1$DE == "Up"]),
  cbind(set = "YBR2_BY4742_ATG_vs_YBR2_BY4742_WT_DOWN",
        orf_name = dea_results.lfc1$orf_name[dea_results.lfc1$contrast == "YBR2_BY4742_ATG_vs_YBR2_BY4742_WT" &
                                               dea_results.lfc1$DE == "Down"]),
  cbind(set = "YBR2_BY4742_ATG_vs_YBR2_BY4742_WT_UP",
        orf_name = dea_results.lfc1$orf_name[dea_results.lfc1$contrast == "YBR2_BY4742_ATG_vs_YBR2_BY4742_WT" &
                                               dea_results.lfc1$DE == "Up"]),
  cbind(set = "YBR2_BY4742_DEL_vs_YBR2_BY4742_WT_DOWN",
        orf_name = dea_results.lfc1$orf_name[dea_results.lfc1$contrast == "YBR2_BY4742_DEL_vs_YBR2_BY4742_WT" &
                                               dea_results.lfc1$DE == "Down"]),
  cbind(set = "YBR2_BY4742_DEL_vs_YBR2_BY4742_WT_UP",
        orf_name = dea_results.lfc1$orf_name[dea_results.lfc1$contrast == "YBR2_BY4742_DEL_vs_YBR2_BY4742_WT" &
                                               dea_results.lfc1$DE == "Up"]),
  cbind(set = "YBR3_BY4741_DEL_vs_YBR3_BY4741_WT_DOWN",
        orf_name = dea_results.lfc1$orf_name[dea_results.lfc1$contrast == "YBR3_BY4741_DEL_vs_YBR3_BY4741_WT" &
                                               dea_results.lfc1$DE == "Down"]),
  cbind(set = "YBR3_BY4741_DEL_vs_YBR3_BY4741_WT_UP",
        orf_name = dea_results.lfc1$orf_name[dea_results.lfc1$contrast == "YBR3_BY4741_DEL_vs_YBR3_BY4741_WT" &
                                               dea_results.lfc1$DE == "Up"]),
  cbind(set = "YBR3_BY4742_ATG_vs_YBR3_BY4742_WT_DOWN",
        orf_name = dea_results.lfc1$orf_name[dea_results.lfc1$contrast == "YBR3_BY4742_ATG_vs_YBR3_BY4742_WT" &
                                               dea_results.lfc1$DE == "Down"]),
  cbind(set = "YBR3_BY4742_ATG_vs_YBR3_BY4742_WT_UP",
        orf_name = dea_results.lfc1$orf_name[dea_results.lfc1$contrast == "YBR3_BY4742_ATG_vs_YBR3_BY4742_WT" &
                                               dea_results.lfc1$DE == "Up"])
)

dea_results.lfc1.orf_sets %>%
  data.frame() %>% group_by(set) %>% count() %>% data.frame()


# ##### MEASURE OVERLAP BETWEEN VARIOUS REPLICATES FORM THE MASS SPEC DATA
# unique(pddat$ID)
# 
combine_replicates <- 'Y'

test_group <- c("attempt_2_WT_YBR_mNG", "attempt_2_doa10D_YBR_mNG", "attempt_1_doa10D_YBR_mNG", "attempt_2_WT_YBR_HA")
control_group <- c("attempt_2_WT_mNG", "attempt_2_doa10D_mNG", "attempt_1_WT_mNG", "attempt_2_WT_HA")

resres <- ms_overlap(combine_replicates, test_group, control_group)

# upset_plot <- resres[[1]]
# result_summary <- resres[[2]]
# result_proteins <- resres[[3]]
# overlap_genes <- resres[[4]]
result_all_proteins <- resres[[5]]
# results_go <- resres[[6]]
# results_kegg <- resres[[7]]

overlap_gene.descriptions <- all_genes[(all_genes$ORF %in% result_all_proteins$Gene) |
                                         (all_genes$GENENAME %in% result_all_proteins$Gene),]
result_all_proteins <- rbind(merge(result_all_proteins, overlap_gene.descriptions, by.x = 'Gene', by.y = 'GENENAME')[c(1:10,15)],
      merge(result_all_proteins, overlap_gene.descriptions, by.x = 'Gene', by.y = 'ORF')[c(1:10,15)])
result_all_proteins <- result_all_proteins[order(result_all_proteins$lfc, result_all_proteins$peptide_cnts_test, decreasing = T),]
rownames(result_all_proteins) <- NULL
write.csv(result_all_proteins, file = 'output/pull_down_two_mNG_samples_overlap.csv', row.names = F)


##### LOAD ALL OVERLAP RESULTS
overlap_gene.one_two_all <- read.csv(file = 'output/pull_down_one_and_two_all_samples_overlap.csv')
overlap_gene.one_two_mng <- read.csv(file = 'output/pull_down_one_and_two_mNG_samples_overlap.csv')
overlap_gene.two_ha_mng <- read.csv(file = 'output/pull_down_two_HA_and_mNG_samples_overlap.csv')
overlap_gene.two_ha <- read.csv(file = 'output/pull_down_two_HA_samples_overlap.csv')
overlap_gene.two_mng <- read.csv(file = 'output/pull_down_two_mNG_samples_overlap.csv')

overlap_gene.all <- rbind(
  cbind(set = "COIP_ALL", orf_name = all_genes$ORF[(all_genes$ORF %in% overlap_gene.one_two_all$Gene[overlap_gene.one_two_all$lfc > 1 & overlap_gene.one_two_all$test_overlap == 1]) |
                                                     (all_genes$GENENAME %in% overlap_gene.one_two_all$Gene[overlap_gene.one_two_all$lfc > 1 & overlap_gene.one_two_all$test_overlap == 1])]),
  cbind(set = "COIP_mNG", orf_name = all_genes$ORF[(all_genes$ORF %in% overlap_gene.one_two_mng$Gene[overlap_gene.one_two_mng$lfc > 1 & overlap_gene.one_two_mng$test_overlap == 1]) |
                                                     (all_genes$GENENAME %in% overlap_gene.one_two_mng$Gene[overlap_gene.one_two_mng$lfc > 1 & overlap_gene.one_two_mng$test_overlap == 1])]),
  cbind(set = "COIP_TWO_ALL", orf_name = all_genes$ORF[(all_genes$ORF %in% overlap_gene.two_ha_mng$Gene[overlap_gene.two_ha_mng$lfc > 1 & overlap_gene.two_ha_mng$test_overlap == 1]) |
                                                     (all_genes$GENENAME %in% overlap_gene.two_ha_mng$Gene[overlap_gene.two_ha_mng$lfc > 1 & overlap_gene.two_ha_mng$test_overlap == 1])]),
  cbind(set = "COIP_TWO_HA", orf_name = all_genes$ORF[(all_genes$ORF %in% overlap_gene.two_ha$Gene[overlap_gene.two_ha$lfc > 1 & overlap_gene.two_ha$test_overlap == 1]) |
                                                     (all_genes$GENENAME %in% overlap_gene.two_ha$Gene[overlap_gene.two_ha$lfc > 1 & overlap_gene.two_ha$test_overlap == 1])]),
  cbind(set = "COIP_TWO_mNG", orf_name = all_genes$ORF[(all_genes$ORF %in% overlap_gene.two_mng$Gene[overlap_gene.two_mng$lfc > 1 & overlap_gene.two_mng$test_overlap == 1]) |
                                                     (all_genes$GENENAME %in% overlap_gene.two_mng$Gene[overlap_gene.two_mng$lfc > 1 & overlap_gene.two_mng$test_overlap == 1])])
  )


##### COIP OVERLAPPING GENES AND VARIOUS ORFSETS
overlap_gene.all %>%
  data.frame() %>%
  filter(orf_name %in% get_snd_genes)

overlap_gene.all %>%
  data.frame() %>%
  filter(orf_name %in% vma_genes)

ybr_dea_coip_overlap <- merge(overlap_gene.all %>% data.frame(),
      dea_results.lfc1.orf_sets %>% data.frame(),
      by = 'orf_name', suffixes = c('_COIP','_DEA')) %>%
  group_by(set_COIP, set_DEA) %>% data.frame()

ybr_dea_coip_overlap.cnts <- merge(merge(ybr_dea_coip_overlap %>%
                                      group_by(set_COIP, set_DEA) %>%
                                      count(),
                                    overlap_gene.all %>% data.frame() %>%
                                      group_by(set) %>%
                                      count(),
                                    by.x = 'set_COIP', by.y = 'set',
                                    suffixes = c('','_total_COIP')),
                              dea_results.lfc1.orf_sets %>% data.frame() %>%
                                group_by(set) %>%
                                count(),
                              by.x = 'set_DEA', by.y = 'set',
                              suffixes = c('_overlap','_total_DEA')) %>%
  mutate(jaccard = n_overlap/(n_total_COIP+n_total_DEA-n_overlap))
ybr_dea_coip_overlap.cnts$max_overlap <- apply(ybr_dea_coip_overlap.cnts[,c('n_total_COIP','n_total_DEA')], 1, FUN = min, na.rm = TRUE)
ybr_dea_coip_overlap.cnts$overlap_coef <- ybr_dea_coip_overlap.cnts$n_overlap/ybr_dea_coip_overlap.cnts$max_overlap
ybr_dea_coip_overlap.cnts <- ybr_dea_coip_overlap.cnts[order(ybr_dea_coip_overlap.cnts$overlap_coef, decreasing = T),]
ybr_dea_coip_overlap.cnts %>%
  filter(n_overlap > 1)

write.csv(all_genes %>%
  filter(ORF %in% ybr_dea_coip_overlap$orf_name[ybr_dea_coip_overlap$set_COIP == 'COIP_mNG' &
                                                  ybr_dea_coip_overlap$set_DEA == 'YBR2_BY4742_ATG_vs_YBR2_BY4742_WT_DOWN']),
  file = 'output/pull_down_mNG_YBR2_BY4742_ATG_vs_YBR2_BY4742_WT_DOWN_DEA_overlap.csv',
  row.names = F)

(ybr_dea_coip_overlap$orf_name[ybr_dea_coip_overlap$set_COIP == 'COIP_mNG' &
                                ybr_dea_coip_overlap$set_DEA == 'YBR2_BY4742_DEL_vs_YBR2_BY4742_WT_DOWN'] %in%
ybr_dea_coip_overlap$orf_name[ybr_dea_coip_overlap$set_COIP == 'COIP_mNG' &
                                ybr_dea_coip_overlap$set_DEA == 'YBR2_BY4742_ATG_vs_YBR2_BY4742_WT_DOWN'])


##### COEXP RESULTS FOR THE OVERLAP GENES
coexp.ybr %>%
  filter(orf_name %in% overlap_gene.one_two_all$ORF)

coexp.ybr %>%
  filter(orf_name %in% get_snd_genes)

##### GSEA WITH COEXP OF THE OVERLAP RESULTS
coexp.annotated.ybr <- coexp.ybr %>% filter(!is.na(orf_name))
ybr_coexp_universe <- coexp.annotated.ybr$cor[!is.na(coexp.annotated.ybr$ENSEMBL)]
names(ybr_coexp_universe) <- coexp.annotated.ybr$ENSEMBL[!is.na(coexp.annotated.ybr$ENSEMBL)]

gsea_ybr <- GSEA(ybr_coexp_universe,
                 TERM2GENE = rbind(orfsets,
                                   data.frame(rbind(dea_results.lfc1.orf_sets,
                                                    overlap_gene.all,
                                                    cbind(set = 'JPN', orf_name = jap_genes),
                                                    cbind(set = 'VMA', orf_name = vma_genes),
                                                    cbind(set = 'GET/SND', orf_name = get_snd_genes)))),
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

write.csv(gsea_ybr, file = 'output/pull_down_coexp_gsea.csv', row.names = F)


##### OVERLAP GENE GO DETAILS
merge(merge(overlap_gene.all %>%
        data.frame() %>% filter(set == 'COIP_TWO_HA'),
      orfsets, by = 'orf_name', suffixes = c('_COIP','_GO')),
      godetails, by.x = 'set_GO', by.y = 'ID') %>%
  group_by(set_GO, Term) %>% count() %>% filter(n > 10) %>% data.frame()


merge(merge(overlap_gene.all %>%
              data.frame() %>% filter(set == 'COIP_TWO_HA'),
            orfsets, by = 'orf_name', suffixes = c('_COIP','_GO')),
      godetails, by.x = 'set_GO', by.y = 'ID') %>%
  group_by(orf_name) %>% count() %>% filter(n > 10) %>% data.frame()

