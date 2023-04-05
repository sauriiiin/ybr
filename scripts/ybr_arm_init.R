##### HA 10/36 ARM Repeat Proteins
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 02/16/2023

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

source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

##### GATHER DATA
all_genes <- dbGetQuery(conn, 'select * from RNASEQ_YBRATG_ALL_GENES')

reg_path <- "/home/sbp29/R/Projects/ybr/input/ms_data/HA_36/arm_regulators/"
reg_files <- list.files(path = reg_path) 

phe_path <- "/home/sbp29/R/Projects/ybr/input/ms_data/HA_36/arm_phenome/"
phe_files <- list.files(path = phe_path) 

yeast_tf <- read.csv(file = '/home/sbp29/R/Projects/ybr/input/ms_data/TF_consensus.csv')
yeast_tf <- str_to_upper(unique(strtrim(yeast_tf$Transcription.Factor,4)))
# yeast_tf[yeast_tf %notin% unique(all_genes$GENENAME)]
yeast_tf <- all_genes %>% filter(GENENAME %in% yeast_tf)

coexp <- readRDS('/home/sbp29/R/Data/20221122_coexpression.RDS')


##### PROCESS COEXP DATA
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


##### GATHER DATA ON REGULATORS
arm.regs <- NULL
for (f in reg_files) {
 temp <- read.table(sprintf('%s%s',reg_path,f), skip = 8, sep = '\t', header = T) 
 arm.regs <- rbind(arm.regs, temp[,c(1:5)])
}
arm.regs <- data.frame(arm.regs)

arm.regs %>%
  group_by(Regulator, Regulator.Systematic.Name) %>%
  count() %>% data.frame()


##### ARM REGS AND THE COEXPRESSION NETWORK
quantile(coexp.ybr$cor, c(0.5,0.95,0.99))

coexp.ybr %>%
  filter(orf_name %in% unique(arm.regs$Regulator.Systematic.Name))

coexp.ybr %>%
  filter(orf_name %in% unique(arm.regs$Target.Systematic.Name))


head(coexp)
# coexp.ybr.arm.int <- coexp %>%
#   filter((gene1 %in% unique(arm.regs$Regulator.Systematic.Name) |
#             gene1 %in% unique(arm.regs$Target.Systematic.Name)) &
#            (gene2 %in% unique(arm.regs$Regulator.Systematic.Name) |
#               gene2 %in% unique(arm.regs$Target.Systematic.Name)))
# 
# ybr.arm.int <- data.frame(rbind(cbind(node = 'Regulator', Gene = unique(arm.regs$Regulator)),
#                                 cbind(node = 'ARM_interactor', Gene = unique(arm.regs$Target))))
# ybr.arm.int <- merge(ybr.arm.int, all_genes, by.x = 'Gene', by.y = 'GENENAME')

coexp.ybr.arm.int <- coexp %>%
  filter((gene1 %in% yeast_tf$ORF |
            gene1 %in% unique(arm.regs$Target.Systematic.Name)) &
           (gene2 %in% yeast_tf$ORF |
              gene2 %in% unique(arm.regs$Target.Systematic.Name)))

ybr.arm.int <- data.frame(rbind(cbind(node = 'Regulator', Gene = yeast_tf$GENENAME),
                                cbind(node = 'ARM_interactor', Gene = unique(arm.regs$Target))))
ybr.arm.int <- merge(ybr.arm.int, all_genes, by.x = 'Gene', by.y = 'GENENAME')
# write.csv(ybr.arm.int, file = 'output/pull_down_two_ha_edger_arm_proteins_regulators.csv', row.names = F)


coexp.ybr.arm.int <- rbind(coexp.ybr.arm.int, coexp.ybr[coexp.ybr$orf_name %in% ybr.arm.int$ORF, colnames(coexp.ybr.arm.int)])
coexp.ybr.arm.int <- merge(merge(coexp.ybr.arm.int, rbind(ybr.arm.int[,c('ORF','Gene')], c('YBR196C-A','YBR196C-A')), by.x = 'gene1', by.y = 'ORF'),
      rbind(ybr.arm.int[,c('ORF','Gene')], c('YBR196C-A','YBR196C-A')), by.x = 'gene2', by.y = 'ORF', suffixes = c('_1','_2'))
# write.csv(coexp.ybr.arm.int, file = 'output/pull_down_two_ha_edger_arm_proteins_regulators_coexp.csv', row.names = F)

coexp.ybr.arm.int %>% filter(cor > 0.7) %>% nrow()
# write.csv(coexp.ybr.arm.int %>% filter(cor > 0.7), file = 'output/pull_down_two_ha_edger_arm_proteins_regulators_coexp_07.csv', row.names = F)

head(coexp.ybr.arm.int)
coexp.ybr.arm.int %>%
  filter(Gene_1 %in% unique(arm.regs$Target),
         Gene_2 %in% unique(arm.regs$Target))


##### GATHER DATA ON ARM PHENOME
arm.phe <- NULL
for (f in phe_files) {
  temp <- read.table(sprintf('%s%s',phe_path,f), sep = '\t', header = T)
  temp$arm_protein <- str_remove(f, '.txt')
  arm.phe <- rbind(arm.phe, temp)
}
arm.phe <- data.frame(arm.phe)

arm.phe %>%
  group_by(arm_protein) %>%
  summarise(p99 = quantile(Correlation.mean, c(0.99)), .groups = 'keep')

merge(arm.phe, arm.phe %>%
  group_by(arm_protein) %>%
  summarise(p99 = quantile(Correlation.mean, c(0.99)), .groups = 'keep'),
  by = 'arm_protein') %>%
  filter(Correlation.mean > p99, Gene.common.name %in% ybr.arm.int$Gene)

arm.phe[,c('Correlation.mean','Gene.common.name','arm_protein')] %>%
  mutate(arm_protein1 = arm_protein, arm_protein2 = Gene.common.name, correlation = Correlation.mean) %>%
  filter(Gene.common.name %in% ybr.arm.int$Gene[ybr.arm.int$node == 'ARM_interactor'])


###### KAP123 COEXP
coexp.kap123 <- coexp %>%
  filter((gene1 == 'YER110C') |
           (gene2 == 'YER110C'))
coexp.kap123$orf_name[!is.na(coexp.kap123$gene1) & coexp.kap123$gene1 == 'YER110C'] <- coexp.kap123$gene2[!is.na(coexp.kap123$gene1) & coexp.kap123$gene1 == 'YER110C']
coexp.kap123$orf_name[!is.na(coexp.kap123$gene2) & coexp.kap123$gene2 == 'YER110C'] <- coexp.kap123$gene1[!is.na(coexp.kap123$gene2) & coexp.kap123$gene2 == 'YER110C']
coexp.kap123 <- merge(coexp.kap123, all_genes, by.x = 'orf_name', by.y = 'ORF', all.x = T)
coexp.kap123 <- coexp.kap123[order(coexp.kap123$cor, decreasing = T),]
coexp.kap123 <- coexp.kap123 %>% filter(!is.na(orf_name))
row.names(coexp.kap123) <- NULL

coexp.kap123 %>% filter(GENENAME %in% c('IRE1','HAC1'))

coexp.kap123.ybr <- merge(coexp.kap123[,c(-2,-3)] %>% filter(!is.na(orf_name)),
                        coexp.ybr[,c(-2,-3)] %>% filter(!is.na(orf_name)),
                        by = c('orf_name','SGD','ENTREZID','GENENAME','ENSEMBL','DESCRIPTION'),
                        suffixes = c('_kap123','_ybr'))

coexp.kap123.ybr %>%
  ggplot(aes(x = cor_kap123, y = cor_ybr)) +
  geom_point() +
  stat_cor(method = 'spearman')


##### GCN4 COEXP
coexp.gcn4 <- coexp %>%
  filter(gene1 == 'YEL009C' | gene2 == 'YEL009C')
coexp.gcn4$orf_name[!is.na(coexp.gcn4$gene1) & coexp.gcn4$gene1 == 'YEL009C'] <- coexp.gcn4$gene2[!is.na(coexp.gcn4$gene1) & coexp.gcn4$gene1 == 'YEL009C']
coexp.gcn4$orf_name[!is.na(coexp.gcn4$gene2) & coexp.gcn4$gene2 == 'YEL009C'] <- coexp.gcn4$gene1[!is.na(coexp.gcn4$gene2) & coexp.gcn4$gene2 == 'YEL009C']
coexp.gcn4 <- merge(coexp.gcn4, all_genes, by.x = 'orf_name', by.y = 'ORF', all.x = T)
coexp.gcn4 <- coexp.gcn4[order(coexp.gcn4$cor, decreasing = T),]
row.names(coexp.gcn4) <- NULL

head(coexp.gcn4)

coexp.gcn4.ybr <- merge(coexp.gcn4[,c(-2,-3)] %>% filter(!is.na(orf_name)),
      coexp.ybr[,c(-2,-3)] %>% filter(!is.na(orf_name)),
      by = c('orf_name','SGD','ENTREZID','GENENAME','ENSEMBL','DESCRIPTION'),
      suffixes = c('_gcn4','_ybr'))
coexp.gcn4.ybr <- coexp.gcn4.ybr[order(coexp.gcn4.ybr$cor_ybr, decreasing = T),]
row.names(coexp.gcn4.ybr) <- NULL

coexp.gcn4.ybr %>%
  ggplot(aes(x = cor_gcn4, y = cor_ybr)) +
  geom_point() +
  stat_cor(method = 'spearman') +
  coord_cartesian(xlim = c(-0.4,0.9),
                  ylim = c(-0.4,0.9))


##### CORRELATION OF COEXPRESSION YBR VS ALL ANNOTATED ORFS
coexp.ybr.annotated <- coexp.ybr %>% filter(!is.na(orf_name))
row.names(coexp.ybr.annotated) <- NULL

coexp.corr.ybr <- NULL
for (o in coexp.ybr.annotated$orf_name) {
  coexp.temp <- coexp %>% filter(gene1 == o | gene2 == o)
  coexp.temp$orf_name[!is.na(coexp.temp$gene1) & coexp.temp$gene1 == o] <- coexp.temp$gene2[!is.na(coexp.temp$gene1) & coexp.temp$gene1 == o]
  coexp.temp$orf_name[!is.na(coexp.temp$gene2) & coexp.temp$gene2 == o] <- coexp.temp$gene1[!is.na(coexp.temp$gene2) & coexp.temp$gene2 == o]
  coexp.temp <- coexp.temp %>% filter(!is.na(orf_name))

  coexp.ybr.orf <- merge(coexp.temp[,c('orf_name','cor')],
                          coexp.ybr[,c('orf_name','cor')],
                          by = 'orf_name',
                          suffixes = c('_orf','_ybr'))
  
  temp.coexp.corr <- coexp.ybr.annotated$cor[coexp.ybr.annotated$orf_name == o]
  temp.corr <- cor(coexp.ybr.orf$cor_ybr, coexp.ybr.orf$cor_orf, method = 'spearman')

  coexp.corr.ybr <- rbind(coexp.corr.ybr, data.frame(orf_name = o,
                                                     cor = temp.coexp.corr,
                                                     cor_cor = temp.corr))
}

coexp.corr.ybr <- merge(coexp.corr.ybr, all_genes, by.x = 'orf_name', by.y = 'ORF')
coexp.corr.ybr <- coexp.corr.ybr[order(coexp.corr.ybr$cor, decreasing = T),]
row.names(coexp.corr.ybr) <- NULL

head(coexp.corr.ybr)
coexp.corr.ybr.tf <- coexp.corr.ybr %>%
  filter(GENENAME %in% yeast_tf)


coexp.corr.ybr %>%
  ggplot(aes(x = cor, y = cor_cor)) +
  geom_point() +
  stat_cor(method = 'spearman') +
  geom_point(data = coexp.corr.ybr.tf,
             aes(x = cor, y = cor_cor), col = 'red') +
  geom_point(data = coexp.corr.ybr %>% filter(GENENAME %in% ybr_ha_ppi$orf_name),
             aes(x = cor, y = cor_cor), col = 'blue')


##### ALL YEAST PROTEINS WITH ARM DOMAINS
all_arm_proteins <- read.csv(file = 'input/ms_data/SSF48371_Genes.csv')

coexp.ybr.arm.int <- coexp %>%
  filter((gene1 %in% all_arm_proteins$Gene.Systematic.Name) &
           (gene2 %in% all_arm_proteins$Gene.Systematic.Name))

coexp.ybr.arm.int <- rbind(coexp.ybr.arm.int, coexp.ybr[coexp.ybr$orf_name %in% all_arm_proteins$Gene.Systematic.Name, colnames(coexp.ybr.arm.int)])

coexp.ybr.arm <- coexp.ybr.arm.int %>% filter(cor > 0.9)
coexp.ybr.arm$row <- coexp.ybr.arm$gene1
coexp.ybr.arm$column <- coexp.ybr.arm$gene2
coexp.ybr.arm <- merge(merge(coexp.ybr.arm, all_arm_proteins[,c(1,2)], by.x = 'row', by.y = 'Gene.Systematic.Name'),
      all_arm_proteins[,c(1,2)], by.x = 'column', by.y = 'Gene.Systematic.Name', suffixes = c('_1','_2'))
coexp.ybr.arm$gene1 <- coexp.ybr.arm$Gene.Name_1
coexp.ybr.arm$gene2 <- coexp.ybr.arm$Gene.Name_2
coexp.ybr.arm <- coexp.ybr.arm[,c(-7,-8)]
write.csv(coexp.ybr.arm, file = 'output/all_arm_proteins_coexp.csv')

all_arm_proteins$ybr_interactor[all_arm_proteins$Gene.Name %in% hello$orf_name] <- 'Yes'
all_arm_proteins$ybr_interactor[is.na(all_arm_proteins$ybr_interactor)] <- 'No'

all_arm_proteins %>% filter(ybr_interactor == 'Yes') %>% count()
write.csv(all_arm_proteins, file = 'output/all_arm_proteins.csv')


paste(all_arm_proteins$Gene.Systematic.Name[all_arm_proteins$ybr_interactor == 'Yes' & all_arm_proteins$Gene.Name %notin% c('LOS1','MTR10','RPN2')],
      collapse = ',')
