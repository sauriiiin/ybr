source('/home/sbp29/R/Projects/methionine/paper/scripts/initialize.R')
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")
orfid_crossmap <- dbGetQuery(conn, 'select * from CARVUNIS_YEAST.orf_ids_crossmap')
orfid_crossmap <- orfid_crossmap[,c(1,8)]
colnames(orfid_crossmap) <- c('gene_id', 'orf_id')
head(orfid_crossmap)

orfs <- readRDS('/home/aar75/rna_seq/Salmon_4_6_22/translatome.rds')

##### FEATURE COUNTS
library(Rsubread)
annot <- read.table(file = '/home/sbp29/RNASeq/YBR/translated_orfs_with_transient_status.gff3')
annot <- annot[,c(9,1,4,5,7)]
annot$V9 <- str_remove(str_split(annot$V9, ';', simplify = T)[,1], 'ID=')
colnames(annot) <- c('GeneID','Chr','Start','End','Strand')
annot$Chr[annot$Chr == 'chrI'] <- 'chr1'
annot$Chr[annot$Chr == 'chrII'] <- 'chr2'
annot$Chr[annot$Chr == 'chrIII'] <- 'chr3'
annot$Chr[annot$Chr == 'chrIV'] <- 'chr4'
annot$Chr[annot$Chr == 'chrV'] <- 'chr5'
annot$Chr[annot$Chr == 'chrVI'] <- 'chr6'
annot$Chr[annot$Chr == 'chrVII'] <- 'chr7'
annot$Chr[annot$Chr == 'chrVIII'] <- 'chr8'
annot$Chr[annot$Chr == 'chrIX'] <- 'chr9'
annot$Chr[annot$Chr == 'chrX'] <- 'chr10'
annot$Chr[annot$Chr == 'chrXI'] <- 'chr11'
annot$Chr[annot$Chr == 'chrXII'] <- 'chr12'
annot$Chr[annot$Chr == 'chrXIII'] <- 'chr13'
annot$Chr[annot$Chr == 'chrXIV'] <- 'chr14'
annot$Chr[annot$Chr == 'chrXV'] <- 'chr15'
annot$Chr[annot$Chr == 'chrXVI'] <- 'chr16'

sample_30b <- read.delim('/home/aar75/rna_seq/longreads/nanopore_ybr/20220819/30b.bed',header=F,col.names=c('chr','start','stop','name','score','strand'))
head(sample_30b)
unique(sample_30b$chr)

sample_30b %>%
  group_by(chr) %>%
  count()

fc <- featureCounts(files=c('/home/aar75/rna_seq/longreads/nanopore_ybr/20220819/30a.sam',
                               '/home/aar75/rna_seq/longreads/nanopore_ybr/20220819/30b.sam',
                            '/home/aar75/rna_seq/longreads/nanopore_ybr/20220819/30c.sam',
                            '/home/aar75/rna_seq/longreads/nanopore_ybr/20220819/30d.sam'),
                    annot.ext=annot,
                    allowMultiOverlap = TRUE,
                    # isPairedEnd=FALSE,
                    # countMultiMappingReads=TRUE,
                    isLongRead = TRUE,
                    nthreads=12)

head(fc$counts)
fc$counts %>%
  melt() %>%
  group_by(Var2) %>%
  summarize(s = sum(value))


reads <- fc$counts %>% data.frame()
colnames(reads) <- c('crispy','wt','stop','syn')
reads$gene_id <- rownames(reads)
rownames(reads) <- NULL

reads <- merge(reads, orfid_crossmap, by = 'gene_id')
reads <- merge(orfs, reads, by = 'orf_id', all.y = T)

reads[,c('orf_id','chr','crispy','wt','stop','syn')] %>%
  melt(id.vars = c('orf_id','chr')) %>%
  filter(variable == 'wt') %>%
  group_by(variable, chr) %>%
  summarize(s = sum(value)) %>%
  data.frame()

reads %>%
  filter(gene_systematic_name == 'YBR196C-A')

reads$crispy[reads$gene_systematic_name == 'YBR196C-A' & !is.na(reads$gene_systematic_name)] <- 0

##### EDGER
library(edgeR)
cnts <- reads[!is.na(reads$transcript), c('wt','crispy','syn','stop')]
rownames(cnts) <- reads$transcript[!is.na(reads$transcript)]

cnts <- cnts[rowSums(cnts) > 1,]

### NORMALIZED COUNTS
dds <- DGEList(cnts, group = rep(1:4,each=1))
dds <- calcNormFactors(dds, method='TMM')
plotMDS(dds)

cnts_tmm <- cpm(dds)
reads_nrm <- cnts_tmm %>% data.frame()
reads_nrm$transcript <- rownames(reads_nrm)
rownames(reads_nrm) <- NULL
head(reads_nrm)

reads_nrm <- merge(orfs[,c('transcript','gene_systematic_name')], reads_nrm, by = 'transcript')

reads_nrm %>%
  filter(gene_systematic_name == 'YBR196C-A')

# reads_nrm <- reads_nrm %>%
#   mutate(crispy_diff = crispy/wt, stop_diff = stop/syn, syn_diff = syn/wt)

### DEG
res <- NULL

temp <- exactTest(dds, pair = c(1,2), dispersion = 0.2^2)
reads_nrm$crispy_vs_wt_lfc <- temp$table$logFC
reads_nrm$crispy_vs_wt_p <- temp$table$PValue
temp$table$padj <- p.adjust(temp$table$PValue, method="BH")
temp$table$transcript <- rownames(temp$table)
rownames(temp$table) <- NULL
res <- rbind(res, data.frame(temp$table, comp = 'CRISPY_vs_WT'))

temp <- exactTest(dds, pair = c(3,4), dispersion = 0.2^2) 
reads_nrm$stop_vs_syn_lfc <- temp$table$logFC
reads_nrm$stop_vs_syn_p <- temp$table$PValue
temp$table$padj <- p.adjust(temp$table$PValue, method="BH")
temp$table$transcript <- rownames(temp$table)
rownames(temp$table) <- NULL
res <- rbind(res, data.frame(temp$table, comp = 'STOP_vs_SYN'))

temp <- exactTest(dds, pair = c(1,4), dispersion = 0.2^2)
reads_nrm$stop_vs_wt_lfc <- temp$table$logFC
reads_nrm$stop_vs_wt_p <- temp$table$PValue
temp$table$padj <- p.adjust(temp$table$PValue, method="BH")
temp$table$transcript <- rownames(temp$table)
rownames(temp$table) <- NULL
res <- rbind(res, data.frame(temp$table, comp = 'STOP_vs_WT'))

temp <- exactTest(dds, pair = c(1,3), dispersion = 0.2^2)
reads_nrm$syn_vs_wt_lfc <- temp$table$logFC
reads_nrm$syn_vs_wt_p <- temp$table$PValue
temp$table$padj <- p.adjust(temp$table$PValue, method="BH")
temp$table$transcript <- rownames(temp$table)
rownames(temp$table) <- NULL
res <- rbind(res, data.frame(temp$table, comp = 'SYN_vs_WT'))

head(res)
merge(res, orfs, by = 'transcript') %>%
  filter(gene_systematic_name == 'YBR196C-A')

hi <- merge(res, orfs, by = 'transcript') %>%
  filter(logFC < 1, PValue <= 0.05, gene_systematic_name != 'X', comp == 'CRISPY_vs_WT')
paste0(hi$gene_systematic_name, collapse = ', ')

res %>%
  ggplot(aes(x = logFC, y = -log(PValue,10))) +
  geom_point() +
  geom_hline(yintercept = -log(0.05,10)) +
  # scale_y_log10() +
  facet_wrap(.~comp)

reads_nrm <- merge(reads_nrm, reads[,c('transcript','crispy','wt','stop','syn')],
                   by = 'transcript', suffixes = c('','_raw'))
head(reads_nrm)

reads_nrm %>%
  filter(crispy_vs_wt_p <= 0.05) %>%
  ggplot(aes(x = log(wt,10), y = log(crispy,10))) +
  geom_abline() +
  geom_point(aes(col = crispy_vs_wt_lfc), size = 1) +
  coord_cartesian(xlim = c(2,3),
                  ylim = c(2,3))

reads_nrm %>%
  ggplot(aes(x = log(wt,10), y = log(crispy,10))) +
  geom_point(aes(col = log(crispy/wt,2)))

reads_nrm %>%
  ggplot(aes(x = log(crispy/wt,2), y = log(stop/syn,2))) +
  geom_point()


#####
head(reads_nrm)
reads_nrm %>%
  ggplot(aes(x = seq(1,nrow(reads_nrm)), y = crispy_vs_wt_lfc)) +
  geom_point()

reads_nrm2 <- reads_nrm[,str_detect(colnames(reads_nrm), 'transcript') | str_detect(colnames(reads_nrm), 'gene_systematic_name') | str_detect(colnames(reads_nrm), '_lfc')] %>%
  melt(id.vars = c('transcript','gene_systematic_name'), variable.name = 'comparison', value.name = 'lfc') %>%
  mutate(comparison = str_remove(comparison, '_lfc'))

reads_nrm2 %>%
  group_by(comparison) %>%
  summarize(q = quantile(lfc, 0.975))

reads_nrm2 %>%
  filter(abs(lfc) >= 4, gene_systematic_name != 'X')

cnts %>%
  mutate(transcript = rownames(cnts)) %>%
  filter(transcript %in% c('chr2_612236' ,'chr2_614024' ,'chr2_614215', 'chr2_614702', 'chr2_614936'))

colSums(cnts)

reads_nrm %>%
  filter(transcript %in% c('chr2_612236' ,'chr2_614024' ,'chr2_614215', 'chr2_614702', 'chr2_614936'))

reads_nrm2 %>%
  filter(transcript %in% c('chr2_612236' ,'chr2_614024' ,'chr2_614215', 'chr2_614702', 'chr2_614936'))


orfs %>%
  filter(transcript %in% c('chr2_612236' ,'chr2_614024' ,'chr2_614215', 'chr2_614702', 'chr2_614936'))

orfid_crossmap %>%
  filter(orf_id %in% c(143090,143103,143104,143113,143114))
