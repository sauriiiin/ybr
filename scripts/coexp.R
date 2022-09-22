
source('/home/sbp29/R/Projects/ybr/scripts/initialize.R')

ybr_coip <- read.csv('output/ybr_pd_scores.csv')
colnames(ybr_coip)
ybr_coip <- ybr_coip[,c(-1,-7)]
ybr_coip$coip_rank <- as.numeric(rownames(ybr_coip))

ybr_coex <- readRDS('/home/aar75/coexpression/20220607_ybr_coexpression.RDS')
temp <- bitr(ybr_coex$gene, fromType = "ENSEMBL",
             toType = c("ENTREZID","GENENAME","DESCRIPTION"),
             OrgDb = org.Sc.sgd.db)

ybr_coex2 <- merge(ybr_coex, temp, by.x = 'gene', by.y = 'ENSEMBL', all.x = T)
ybr_coex2 <- ybr_coex2[order(ybr_coex2$rho, decreasing = T),]
rownames(ybr_coex2) <- NULL

head(ybr_coex2)

ybr_coex3 <- ybr_coex2 %>% filter(!is.na(DESCRIPTION))
ybr_coex3$coex_rank <- as.numeric(rownames(ybr_coex3))
ybr_coex3 <- merge(ybr_coex3, ybr_coip, by.x = c('gene', 'GENENAME'), by.y = c('ORF', 'gene'), all.x = T)

ybr_coex3 <- ybr_coex3[order(ybr_coex3$coex_rank, decreasing = F),]
rownames(ybr_coex3) <- NULL

ybr_coex3$ER_related <- str_detect(ybr_coex3$DESCRIPTION, 'ER') | str_detect(ybr_coex3$DESCRIPTION, 'endoplasmic')
ybr_coex3$mito_related <- str_detect(ybr_coex3$DESCRIPTION, 'mitochondria')
ybr_coex3$vac_related <- str_detect(ybr_coex3$DESCRIPTION, 'vacuole')
ybr_coex3$transport_related <- str_detect(ybr_coex3$DESCRIPTION, 'transport')

write.csv(ybr_coex3, file = 'output/ybr_coex_coip.csv', row.names = F)


