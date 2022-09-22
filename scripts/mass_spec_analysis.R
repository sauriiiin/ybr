##### MASS SPEC RESULT ANALYSIS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 08/05/2021 

##### INITIALIZE
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)
library(egg)
library(zoo)
library(ggrepel)
library(ggforce)
library(plotly)
library(scales)
library(reshape2)
library(locfit)
library(growthcurver)
library(rstatix)
library(gtools)
library(locfit)
library(growthrates)
library(RMariaDB)
library(genefilter)
library(apeglm)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(readr)
library(readxl)

out_path <- "~/R/Projects/ybr/output"
fig_path <- "~/R/Projects/ybr/output/figs"

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
source("~/R/Projects/rnaseek/functions/gkEnrich.R")
conn <- initialize.sql("saurin_test")

GFP_crampone <- read_excel("input/interactions/crapome.xlsx")

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

##### LOADING DATA
pro.prbs <- read_csv("input/ms_data/overview_protein_probabilities.csv") %>% data.frame()
temp <- str_split(pro.prbs$Identified.Proteins..219.221., '=', simplify = T)[,3]
temp <- str_remove(temp, ' PE')
pro.prbs$gene <- temp
pro.prbs <- pro.prbs[,c("gene","Molecular.Weight",colnames(pro.prbs)[str_detect(colnames(pro.prbs),'AF')])]
for (i in seq(1,dim(pro.prbs)[[2]])) {
  pro.prbs[,i] <- str_remove(pro.prbs[,i], '%')
}

pep.cnts <- read_csv("input/ms_data/overview_unique_peptide_counts.csv") %>% data.frame()
temp <- str_split(pep.cnts$Identified.Proteins..219.221., '=', simplify = T)[,3]
temp <- str_remove(temp, ' PE')
pep.cnts$gene <- temp
pep.cnts <- pep.cnts[,c("gene","Molecular.Weight",colnames(pep.cnts)[str_detect(colnames(pep.cnts),'AF')])]

ms.dat <- merge(pep.cnts, pro.prbs, by = c('gene','Molecular.Weight'), suffixes = c('_cnt','_prb'))
ms.dat$AFO01_prb <- as.numeric(ms.dat$AFO01_prb)
ms.dat$AFO02_prb <- as.numeric(ms.dat$AFO02_prb)
ms.dat$AFO03_prb <- as.numeric(ms.dat$AFO03_prb)
ms.dat$AFO04_prb <- as.numeric(ms.dat$AFO04_prb)
ms.dat$AFO05_prb <- as.numeric(ms.dat$AFO05_prb)
ms.dat$AFO06_prb <- as.numeric(ms.dat$AFO06_prb)
ms.dat$AFO07_prb <- as.numeric(ms.dat$AFO07_prb)
ms.dat$AFO08_prb <- as.numeric(ms.dat$AFO08_prb)
ms.dat$AFO09_prb <- as.numeric(ms.dat$AFO09_prb)
ms.dat$AFO10_prb <- as.numeric(ms.dat$AFO10_prb)
ms.dat$AFO11_prb <- as.numeric(ms.dat$AFO11_prb)
ms.dat$AFO12_prb <- as.numeric(ms.dat$AFO12_prb)
ms.dat$AFO13_prb <- as.numeric(ms.dat$AFO13_prb)

# spe.cnts <- read_csv("input/ms_data/overview_spectrum_counts.csv") %>% data.frame()
# temp <- str_split(spe.cnts$Identified.Proteins..219.221., '=', simplify = T)[,3]
# temp <- str_remove(temp, ' PE')
# spe.cnts$gene <- temp
# spe.cnts <- spe.cnts[,c("gene","Molecular.Weight",colnames(spe.cnts)[str_detect(colnames(spe.cnts),'AF')])]
# 
# spe.u.cnts <- read_csv("input/ms_data/overview_unique_spectrum_counts.csv") %>% data.frame()
# temp <- str_split(spe.u.cnts$Identified.Proteins..219.221., '=', simplify = T)[,3]
# temp <- str_remove(temp, ' PE')
# spe.u.cnts$gene <- temp
# spe.u.cnts <- spe.u.cnts[,c("gene","Molecular.Weight",colnames(spe.u.cnts)[str_detect(colnames(spe.u.cnts),'AF')])]

##### GENENAME TO ORF
head(ms.dat)

ms.dat <- merge(ms.dat, bitr(ms.dat$gene, fromType = "GENENAME",
                                 toType = c("ORF","DESCRIPTION"),
                                 OrgDb = org.Sc.sgd.db), by.x = 'gene', by.y = 'GENENAME', all = T)
ms.dat$ORF[is.na(ms.dat$ORF)] <- ms.dat$gene[is.na(ms.dat$ORF)]

ms.dat$ybr_cnts <- rowSums(ms.dat[,c(3:8)])
ms.dat$cnt_cnts <- rowSums(ms.dat[,c(9:14)])
ms.dat$cnt_cnts2 <- rowSums(ms.dat[,c(9:15)])

# write.csv(ms.dat, sprintf('%s/ms_dat.csv', out_path))
##### CRAPOME
plot(density(GFP_crampone$NUM_EXPT))
plot(density(rowSums(pep.cnts[3])))

##### BASED ON RAW COUNTS
data.bnd <- data.frame()
data.nodes <- data.frame()
data.nodes.cnts <- data.frame()
goe <- data.frame()
kegg <- data.frame()
for (min_crap in seq(1,15,2)) {
  allgenes <- NULL
  allgenes$ENSEMBL <- keys(org.Sc.sgd.db, column = 'ENSEMBL')
  allgenes$ENSEMBL <- allgenes$ENSEMBL[!(allgenes$ENSEMBL %in% GFP_crampone$GENE[GFP_crampone$NUM_EXPT >= min_crap])]
  for (min_prob in c(80,90)) {
    for (min_cnt in c(1,5)) {
      temp <- ms.dat$gene[((ms.dat$AFO01_cnt >= min_cnt & ms.dat$AFO01_prb > min_prob) |
                             (ms.dat$AFO02_cnt >= min_cnt & ms.dat$AFO02_prb > min_prob) |
                             (ms.dat$AFO03_cnt >= min_cnt & ms.dat$AFO03_prb > min_prob) |
                             (ms.dat$AFO04_cnt >= min_cnt & ms.dat$AFO04_prb > min_prob) |
                             (ms.dat$AFO05_cnt >= min_cnt & ms.dat$AFO05_prb > min_prob) |
                             (ms.dat$AFO06_cnt >= min_cnt & ms.dat$AFO06_prb > min_prob)) &
                            rowSums(ms.dat[9:15]) == 0 &
                            !(ms.dat$gene %in% GFP_crampone$GENE[GFP_crampone$NUM_EXPT >= min_crap])]
      data.bnd <- rbind(data.bnd, data.frame(common = ms.dat$gene, orf_name = ms.dat$ORF,
                                 AFO01 = ms.dat$AFO01_cnt >= min_cnt & ms.dat$AFO01_prb > min_prob,
                                 AFO02 = ms.dat$AFO02_cnt >= min_cnt & ms.dat$AFO02_prb > min_prob,
                                 AFO03 = ms.dat$AFO03_cnt >= min_cnt & ms.dat$AFO03_prb > min_prob,
                                 AFO04 = ms.dat$AFO04_cnt >= min_cnt & ms.dat$AFO04_prb > min_prob,
                                 AFO05 = ms.dat$AFO05_cnt >= min_cnt & ms.dat$AFO05_prb > min_prob,
                                 AFO06 = ms.dat$AFO06_cnt >= min_cnt & ms.dat$AFO06_prb > min_prob,
                                 min_peps = min_cnt, min_prob = min_prob, min_crap_expt = min_crap))
      data.nodes <- rbind(data.nodes, cbind(temp, min_peps = min_cnt, min_prob = min_prob, min_crap_expt = min_crap))
      data.nodes.cnts <- rbind(data.nodes.cnts, cbind(min_peps = min_cnt, min_prob = min_prob, min_crap_expt = min_crap,
                                                just_ybr = length(ms.dat$gene[(ms.dat$AFO01_cnt >= min_cnt & ms.dat$AFO01_prb > min_prob) |
                                                                                (ms.dat$AFO02_cnt >= min_cnt & ms.dat$AFO02_prb > min_prob) |
                                                                                (ms.dat$AFO03_cnt >= min_cnt & ms.dat$AFO03_prb > min_prob) |
                                                                                (ms.dat$AFO04_cnt >= min_cnt & ms.dat$AFO04_prb > min_prob) |
                                                                                (ms.dat$AFO05_cnt >= min_cnt & ms.dat$AFO05_prb > min_prob) |
                                                                                (ms.dat$AFO06_cnt >= min_cnt & ms.dat$AFO06_prb > min_prob)]),
                                                minus_cntrl = length(ms.dat$gene[((ms.dat$AFO01_cnt >= min_cnt & ms.dat$AFO01_prb > min_prob) |
                                                                                    (ms.dat$AFO02_cnt >= min_cnt & ms.dat$AFO02_prb > min_prob) |
                                                                                    (ms.dat$AFO03_cnt >= min_cnt & ms.dat$AFO03_prb > min_prob) |
                                                                                    (ms.dat$AFO04_cnt >= min_cnt & ms.dat$AFO04_prb > min_prob) |
                                                                                    (ms.dat$AFO05_cnt >= min_cnt & ms.dat$AFO05_prb > min_prob) |
                                                                                    (ms.dat$AFO06_cnt >= min_cnt & ms.dat$AFO06_prb > min_prob)) &
                                                                                   rowSums(ms.dat[9:15]) == 0]),
                                                minus_crap = length(ms.dat$gene[((ms.dat$AFO01_cnt >= min_cnt & ms.dat$AFO01_prb > min_prob) |
                                                                                   (ms.dat$AFO02_cnt >= min_cnt & ms.dat$AFO02_prb > min_prob) |
                                                                                   (ms.dat$AFO03_cnt >= min_cnt & ms.dat$AFO03_prb > min_prob) |
                                                                                   (ms.dat$AFO04_cnt >= min_cnt & ms.dat$AFO04_prb > min_prob) |
                                                                                   (ms.dat$AFO05_cnt >= min_cnt & ms.dat$AFO05_prb > min_prob) |
                                                                                   (ms.dat$AFO06_cnt >= min_cnt & ms.dat$AFO06_prb > min_prob)) &
                                                                                  rowSums(ms.dat[9:15]) == 0 &
                                                                                  !(ms.dat$gene %in% GFP_crampone$GENE[GFP_crampone$NUM_EXPT >= min_crap])])))
      
      
      temp <- bitr(temp, fromType = "GENENAME",
                   toType = c("ENTREZID","GENENAME","ENSEMBL"),
                   OrgDb = org.Sc.sgd.db)
      
      temp.goe <- enrichGO(gene          = temp$ENSEMBL,
                           universe      = allgenes$ENSEMBL,
                           OrgDb         = org.Sc.sgd.db,
                           keyType       = "ENSEMBL",
                           ont           = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05)
      if (dim(temp.goe)[1] != 0) {
        goe <- rbind(goe, data.frame(temp.goe, min_peps = min_cnt, min_prob = min_prob, min_crap_expt = min_crap))
      }
      
      temp.kegg <- enrichKEGG(gene         = temp$ENSEMBL,
                              universe     = allgenes$ENSEMBL,
                              organism     = 'sce',
                              pvalueCutoff = 0.05)
      if (dim(temp.kegg)[1] != 0) {
        kegg <- rbind(kegg, data.frame(temp.kegg, min_peps = min_cnt, min_prob = min_prob, min_crap_expt = min_crap))
      }
    }
  }
}
goe$GeneRatio <- as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,2])
goe$BgRatio <- as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,2])
goe$GO <- paste0(goe$ONTOLOGY, '_', goe$Description)

kegg$GeneRatio <- as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,2])
kegg$BgRatio <- as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,2])

goe <- goe[order(goe$min_peps,goe$min_prob,goe$min_crap_expt,-goe$GeneRatio,-goe$Count,goe$qvalue),]
kegg <- kegg[order(kegg$min_peps,kegg$min_prob,kegg$min_crap_expt,-kegg$GeneRatio,-kegg$Count,kegg$qvalue),]

# write.csv(goe, file = sprintf('%s/go_enrichments.csv', out_path), row.names = F)
# write.csv(kegg, file = sprintf('%s/kegg_enrichments.csv', out_path), row.names = F)
write.csv(data.nodes.cnts[order(data.nodes.cnts$min_peps, -data.nodes.cnts$min_prob, data.nodes.cnts$min_crap_expt),],
          file = sprintf('%s/ybr_pi_counts.csv', out_path), row.names = F)

data.nodes <- merge(bitr(data.nodes$temp, fromType = "GENENAME",
     toType = c("ENTREZID","GENENAME","ENSEMBL","DESCRIPTION"),
     OrgDb = org.Sc.sgd.db), data.nodes, by.x = 'GENENAME', by.y = 'temp', all = T)
data.nodes$orf_name <- data.nodes$ENSEMBL
data.nodes$orf_name[is.na(data.nodes$orf_name)] <- data.nodes$GENENAME[is.na(data.nodes$orf_name)]
data.nodes$min_peps <- as.numeric(as.character(data.nodes$min_peps))
data.nodes$min_prob <- as.numeric(as.character(data.nodes$min_prob))
data.nodes$min_crap_expt <- as.numeric(as.character(data.nodes$min_crap_expt))

data.bnd <- data.bnd[data.bnd$min_crap_expt == 1, -11]
data.bnd <- melt(data.bnd, id.vars = c('common','orf_name','min_peps','min_prob'), variable.name = 'band')
data.bnd <- data.bnd[data.bnd$value == T,]
write.csv(data.bnd[,-6], file = sprintf('%s/ms_band_proteins.csv', out_path), row.names = F)

plt.bnds <- data.bnd %>%
  filter(min_peps == 1, min_prob == 90) %>%
  ggplot(aes(x = band, fill = 'Total')) +
  geom_bar(stat = 'count') +
  geom_bar(data = data.bnd %>%
             group_by(common, orf_name, min_peps, min_prob) %>%
             summarize(band = paste(band, collapse = ', '), .groups = 'keep') %>%
             data.frame() %>% filter(nchar(band) <= 5, min_peps == 1, min_prob == 90),
           aes(x = band, fill = 'Unique'), stat = 'count') +
  geom_bar(data = data.bnd[!(data.bnd$orf_name %in% ms.dat$ORF[rowSums(ms.dat[9:15]) != 0]),] %>%
             group_by(common, orf_name, min_peps, min_prob) %>%
             summarize(band = paste(band, collapse = ', '), .groups = 'keep') %>%
             data.frame() %>% filter(nchar(band) <= 5, min_peps == 1, min_prob == 90),
           aes(x = band, fill = 'Unique to YBR'), stat = 'count') +
  scale_fill_manual(name = 'Proteins Per Band',
                      values = c('Total' = '#C5CAE9',
                                 'Unique' = '#448AFF',
                                 'Unique to YBR' = '#303F9F')) +
  labs(title = 'MS Data',
       subtitle = 'Min. No. of Peptides = 1 | Min. Probability = 90%',
       x = 'Gel Band',
       y = 'Count') +
  # facet_wrap(.~min_peps*min_prob) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        plot.subtitle = element_text(size = txt, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))+
  coord_flip() +
  scale_x_discrete(limits = c('AFO06','AFO05','AFO04','AFO03','AFO02','AFO01'))
ggsave(sprintf('%s/MS_BAND_PROTEINS.png',fig_path), plt.bnds,
       width = one.c, height = one.c, units = 'mm',
       dpi = 600)
  

##### INTERACTOR NETWORK
cc.int <- read_excel("input/interactions/cocomplex_interactions.xlsx", col_names = F)
cc.int$type <- 'Co-complex'
y2h.int <- read_excel("input/interactions/Y2H_interactions.xlsx", col_names = F)
y2h.int$type <- 'Y2H'
lc.int <- read_excel("input/interactions/LC_interactions.xlsx", col_names = F)
lc.int$type <- 'Literature'
pre.int <- read_excel("input/interactions/PIPE_interactions.xlsx", col_names = F)
pre.int$type <- 'Predicted'

data.int <- data.frame(rbind(cc.int, y2h.int, lc.int, pre.int))
colnames(data.int) <- c('source','target','interactor')
head(data.int)

pma1.int <- data.int[data.int$source == 'YGL008C',]
data.edges <- NULL
for (min_peps in unique(data.nodes$min_peps)) {
  for (min_prob in unique(data.nodes$min_prob[data.nodes$min_peps == min_peps])) {
    for (min_crap in unique(data.nodes$min_crap_expt[data.nodes$min_prob == min_prob & data.nodes$min_peps == min_peps])) {
      temp <- data.nodes$orf_name[data.nodes$min_peps == min_peps &
                                 data.nodes$min_prob == min_prob &
                                 data.nodes$min_crap_expt == min_crap]
      data.edges <- data.frame(rbind(data.edges, cbind(data.int[data.int$source %in% temp & data.int$target %in% temp,],
                                                 min_peps = min_peps, min_prob = min_prob, min_crap_expt = min_crap)))
      data.nodes$status[data.nodes$min_peps == min_peps &
                       data.nodes$min_prob == min_prob &
                       data.nodes$min_crap_expt == min_crap &
                       data.nodes$orf_name %in% temp[!(temp %in% data.int$target) & !(temp %in% data.int$target)]] <- 'No Interaction Data'
    }
  }
}
data.nodes$status[is.na(data.nodes$status)] <- 'Interaction Data Present'
data.nodes$pma_int[data.nodes$orf_name %in% pma1.int$target] <- 'PMA1 Interactor'
data.nodes$pma_int[is.na(data.nodes$pma_int)] <- 'Not A PMA1 Interactor'
data.nodes <- data.nodes[order(data.nodes$min_peps, data.nodes$min_prob, data.nodes$min_crap_expt),] 
data.nodes <- merge(data.nodes, data.bnd[,c(-1,-6)] %>%
                      group_by(orf_name, min_peps, min_prob) %>%
                      summarize(band = paste(band, collapse = ', '), .groups = 'keep') %>%
                      data.frame() %>% filter(nchar(band) <= 5), 
                    by = c('orf_name','min_peps','min_prob'), all.x = T)

data.edges <- data.edges[data.edges$source != data.edges$target,]
data.edges$pma_int[data.edges$source %in% pma1.int$target] <- 'PMA1 Interactor'
data.edges$pma_int[is.na(data.edges$pma_int)] <- 'Not A PMA1 Interactor'

temp <- bitr(unique(c(data.edges$source, data.edges$target)), fromType = "ENSEMBL",
     toType = "GENENAME", OrgDb = org.Sc.sgd.db, drop = F)
temp$GENENAME[is.na(temp$GENENAME)] <- temp$ENSEMBL[is.na(temp$GENENAME)]
temp <- temp %>% group_by(ENSEMBL) %>%
  summarize(GENENAME = paste(GENENAME, collapse = ',')) %>% data.frame()
temp$GENENAME[temp$ENSEMBL == 'YKR059W'] <- 'TIF1'
temp$GENENAME[temp$ENSEMBL == 'YJL138C'] <- 'TIF2'

data.edges <- merge(temp, data.edges, by.x = 'ENSEMBL', by.y = 'source')
data.edges <- merge(temp, data.edges, by.x = 'ENSEMBL', by.y = 'target')
colnames(data.edges) <- c('target_orf','target','source_orf','source',colnames(data.edges)[5:9])
data.edges <- data.edges[order(data.edges$min_peps, data.edges$min_prob, data.edges$min_crap_expt),]

write.csv(data.nodes, file = sprintf('%s/ybr_nodes.csv', out_path))
write.csv(data.edges, file = sprintf('%s/ybr_edges.csv', out_path))

##### GENEMANIA NETWORK
data.gm <- read.csv(file = 'input/genemania-network-2.txt', header = T, sep = '\t', stringsAsFactors = F)
# data.gm.nodes <- NULL
# data.gm.nodes$genes <- unique(c(unique(data.gm$Entity.1), unique(data.gm$Entity.2)))
# data.gm.nodes <- data.frame(data.gm.nodes, stringsAsFactors = F)
# data.gm.nodes[!(data.gm.nodes$genes %in% data.nodes$GENENAME),]
data.gm %>%
  group_by(Entity.1) %>%
  count() %>% data.frame() %>% arrange(-n) %>%
  head(10)
data.gm %>%
  group_by(Entity.2) %>%
  count() %>% data.frame() %>% arrange(-n)


##### HIGH DEGREE NODES
head(data.edges)
data.edges %>%
  filter(min_peps == 1, min_prob == 90, min_crap_expt == 5) %>%
  group_by(target) %>% count() %>% data.frame()


