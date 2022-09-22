### AFTER SUMMER YBR UPDATE

head(data.nodes)

data.nodes2 <- data.nodes %>%
  filter(min_crap_expt == 15, min_peps == 1, min_prob == 90)
colnames(data.nodes2)

data.nodes2 <- data.nodes2[,c('orf_name','band','GENENAME','DESCRIPTION')]
data.nodes2 <- merge(data.nodes2, data.frame(orf_name = ms.dat$ORF, pep_cnts = rowSums(ms.dat[3:8])), by = 'orf_name')
data.nodes2$crapome[data.nodes2$GENENAME %in% GFP_crampone$GENE[GFP_crampone$NUM_EXPT >= 5]] <- 'Present'
data.nodes2$crapome[is.na(data.nodes2$crapome)] <- 'Absent'

write.csv(data.nodes2, file = 'output/ybr_nodes2.csv')
data.gm.nc <- read.csv('input/genemania-network-nocrap.csv', stringsAsFactors = F)

data.gm.nc$Weight
