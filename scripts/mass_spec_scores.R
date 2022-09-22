

pep.cnts <- read_csv("input/ms_data/overview_unique_peptide_counts.csv") %>% data.frame()
temp <- str_split(pep.cnts$Identified.Proteins..219.221., '=', simplify = T)[,3]
temp <- str_remove(temp, ' PE')
pep.cnts$gene <- temp
pep.cnts <- pep.cnts[,c("gene","Molecular.Weight",colnames(pep.cnts)[str_detect(colnames(pep.cnts),'AF')])]
head(pep.cnts)


pep.cov <- read_csv("input/ms_data/overview_protein_coverage.csv") %>% data.frame()
temp <- str_split(pep.cov$Identified.Proteins..219.221., '=', simplify = T)[,3]
temp <- str_remove(temp, ' PE')
pep.cov$gene <- temp
pep.cov <- pep.cov[,c("gene","Molecular.Weight",colnames(pep.cov)[str_detect(colnames(pep.cov),'AF')])]
head(pep.cov)

ybr.pd <- merge(pep.cnts, pep.cov, by = c('gene','Molecular.Weight'), suffixes = c('_cnt','_cov'))
ybr.pd$pep_cnt <- rowSums(ybr.pd[,c(3:8)])
ybr.pd$pep_cov <- rowMeans(ybr.pd[,c(16:21)], na.rm = T)

ybr.pd <- ybr.pd[,c('gene','pep_cnt','pep_cov')] 
ybr.pd <- merge(ybr.pd, bitr(ybr.pd$gene, fromType = "GENENAME",
                             toType = c("ORF","DESCRIPTION"),
                             OrgDb = org.Sc.sgd.db), by.x = 'gene', by.y = 'GENENAME', all = T)
ybr.pd$ORF[is.na(ybr.pd$ORF)] <-ybr.pd$gene[is.na(ybr.pd$ORF)]


ybr.nodes <- read.csv(file = 'output/ybr_nodes2.csv')
head(ybr.nodes)
gm.nc <- read.csv('input/genemania-network-nocrap.csv', stringsAsFactors = F)
gm.deg <- merge(gm.nc %>%
        group_by(Entity.1) %>%
        count() %>%
        data.frame(),
      gm.nc %>%
        group_by(Entity.2) %>%
        count() %>%
        data.frame(), by.x = 'Entity.1', by.y = 'Entity.2', all = T)
gm.deg$degree <- rowSums(gm.deg[,c(2,3)], na.rm = T)


ybr.pd <- merge(ybr.pd, gm.deg[,c('Entity.1','degree')], by.x = 'gene', by.y = 'Entity.1', all = T)

ybr.pd <- ybr.pd[ybr.pd$ORF %in% ybr.nodes$orf_name,]
ybr.pd$pep_cnt_n <- ybr.pd$pep_cnt/max(ybr.pd$pep_cnt)
ybr.pd$pep_cov_n <- ybr.pd$pep_cov/max(ybr.pd$pep_cov)
ybr.pd$degree_n <- ybr.pd$degree/max(ybr.pd$degree, na.rm = T)

ybr.pd$score <- rowSums(ybr.pd[,c(8:10)], na.rm = T)
ybr.pd <- ybr.pd[order(-ybr.pd$score),]

write.csv(ybr.pd, file = 'output/ybr_pd_scores.csv')
