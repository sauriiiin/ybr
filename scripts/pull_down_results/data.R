

library(dplyr)
library(stringr)


##### PULL DOWN #1
pro.prbs <- read.csv("/home/sbp29/R/Projects/ybr/input/ms_data/overview_protein_probabilities.csv")
temp <- str_split(pro.prbs$Identified.Proteins..219.221., '=', simplify = T)[,3]
temp <- str_remove(temp, ' PE')
pro.prbs$gene <- temp
pro.prbs <- pro.prbs[,c("gene","Molecular.Weight",colnames(pro.prbs)[str_detect(colnames(pro.prbs),'AF')])]
for (i in seq(1,dim(pro.prbs)[[2]])) {
  pro.prbs[,i] <- str_remove(pro.prbs[,i], '%')
}

pep.cnts <- read.csv("/home/sbp29/R/Projects/ybr/input/ms_data/overview_unique_peptide_counts.csv")
temp <- str_split(pep.cnts$Identified.Proteins..219.221., '=', simplify = T)[,3]
temp <- str_remove(temp, ' PE')
pep.cnts$gene <- temp
pep.cnts <- pep.cnts[,c("gene","Molecular.Weight",colnames(pep.cnts)[str_detect(colnames(pep.cnts),'AF')])]
colnames(pep.cnts) <- c('Gene', colnames(pep.cnts)[-1])
pep.cnts <- pep.cnts[,-2]

pep.cnts[,c(2:14)] <- data.frame(pro.prbs[,c(3:15)] %>% mutate_if(is.character, as.numeric) > 90) %>%
  mutate_if(is.logical, as.numeric) * data.frame(pep.cnts[,c(2:14)] %>% mutate_if(is.character, as.numeric))
write.csv(pep.cnts, file = 'scripts/pull_down_results/pull_down_one_cnts.csv', row.names = F)

pep.cnts[,c(2:14)] <- data.frame(pep.cnts[,c(2:14)] > 0) %>%
  mutate_if(is.logical, as.numeric)
write.csv(pep.cnts, file = 'scripts/pull_down_results/pull_down_one.csv', row.names = F)



##### PULL DOWN #2
hatag_yeastonly <- read.csv(file = 'input/ms_data/221225_YBRPullDown_HATag_YeastonlyDB.csv')
hatag_yeastonly <- hatag_yeastonly[hatag_yeastonly$Protein.Probability > 0.9,
                                   c('Gene',colnames(hatag_yeastonly)[str_detect(colnames(hatag_yeastonly), 'AFO')])]

mngtag_yeastonly <- read.csv(file = 'input/ms_data/221225_YBRPullDown_mNGTag_YeastonlyDB.csv')
mngtag_yeastonly <- mngtag_yeastonly[mngtag_yeastonly$Protein.Probability > 0.9,
                                     c('Gene',colnames(mngtag_yeastonly)[str_detect(colnames(mngtag_yeastonly), 'AFO')])]


pdtwo <- merge(hatag_yeastonly, mngtag_yeastonly, by = 'Gene', all = T)
pdtwo[is.na(pdtwo)] <- 0
pdtwo <- pdtwo[pdtwo$Gene != '',]
write.csv(pdtwo, file = 'scripts/pull_down_results/pull_down_two_cnts.csv', row.names = F)

pdtwo[,c(2:81)] <- data.frame(pdtwo[,c(2:81)] > 0) %>% mutate_if(is.logical, as.numeric)
write.csv(pdtwo, file = 'scripts/pull_down_results/pull_down_two.csv', row.names = F)


pdtwo %>%
  filter(Gene == 'GCN1')
pdtwo %>%
  melt(id.vars = 'Gene') %>%
  group_by(variable) %>%
  summarize(total_proteins_found = sum(value)) %>%
  data.frame()
