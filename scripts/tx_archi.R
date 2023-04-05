library(dplyr)
library(RMariaDB)
library(stringr)
library(Rsubread)
library(reshape2)
library(ggplot2)
library(cowplot)

source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

out_path <- "~/R/Projects/ybr/output"


## ORF INFORMATION
orfid_crossmap <- dbGetQuery(conn, 'select * from CARVUNIS_YEAST.orf_ids_crossmap')
orfid_crossmap <- orfid_crossmap[,c(1,8)]
colnames(orfid_crossmap) <- c('gene_id', 'orf_id')

orf_info<-readRDS('/home/aar75/rna_seq/Salmon_4_6_22/translatome.rds')

## ANNOTATION FILE
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

## DIRECT RNASEQ DATA
# sample_30a<-read.delim('/home/aar75/rna_seq/longreads/nanopore_ybr/20220819/30a.bed',header=F,col.names=c('chr','start','stop','name','score','strand'))
sample_30a2<-read.delim('/home/sbp29/RNASeq/YBR/minION/2208/sample_30a_/sample30a.bed',header=F,col.names=c('chr','start','stop','name','score','strand'))
# sample_30b<-read.delim('/home/aar75/rna_seq/longreads/nanopore_ybr/20220819/30b.bed',header=F,col.names=c('chr','start','stop','name','score','strand'))
sample_30b2<-read.delim('/home/sbp29/RNASeq/YBR/minION/2208/sample_30b/sample30b.bed',header=F,col.names=c('chr','start','stop','name','score','strand'))
# sample_30c<-read.delim('/home/aar75/rna_seq/longreads/nanopore_ybr/20220819/30c.bed',header=F,col.names=c('chr','start','stop','name','score','strand'))
sample_30c2<-read.delim('/home/sbp29/RNASeq/YBR/minION/2208/sample_30c/sample30c.bed',header=F,col.names=c('chr','start','stop','name','score','strand'))
# sample_30d<-read.delim('/home/aar75/rna_seq/longreads/nanopore_ybr/20220819/30d.bed',header=F,col.names=c('chr','start','stop','name','score','strand'))
sample_30d2<-read.delim('/home/sbp29/RNASeq/YBR/minION/2208/sample_30d/sample30d.bed',header=F,col.names=c('chr','start','stop','name','score','strand'))

samples <- rbind(
  # cbind(sample_30b, sample = 'wt'),
  cbind(sample_30b2, sample = 'wt2'),
  # cbind(sample_30a, sample = 'crispy'),
  cbind(sample_30a2, sample = 'crispy2'),
  # cbind(sample_30d, sample = 'syn'),
  cbind(sample_30d2, sample = 'syn2'),
  # cbind(sample_30c, sample = 'stop'),
  cbind(sample_30c2, sample = 'stop2')
) %>% data.frame()
# head(samples)

##### COUNTING AND TRANSCRIPT ARCHITECHTURE
loci_range <- seq(1,nrow(annot),1) #seq(14086-1,14086+4,1)
min_overlap <- 10
ybr_transcripts <- samples
ybr_transcripts$orf <- NA
for(i in loci_range){
  temp <- samples[(samples$chr == annot$Chr[i] & samples$strand == annot$Strand[i]) &
                    (((samples$start <= annot$Start[i] & samples$stop >= annot$Start[i]) & (samples$stop - annot$Start[i]) > min_overlap) |
                       ((samples$start <= annot$End[i] & samples$stop >= annot$End[i]) & (annot$End[i] - samples$start) > min_overlap) |
                       ((samples$start >= annot$Start[i] & samples$stop <= annot$End[i]) & (samples$stop - samples$start) > min_overlap) |
                       (samples$start <= annot$Start[i] & samples$stop >= annot$End[i])),]
  nrow(temp)
  if (nrow(temp) > 0) {
    # temp2 <- cbind(temp, hello = 1)
    # colnames(temp2) <- c(colnames(temp2)[-8], annot$GeneID[i])
    # ybr_transcripts <- merge(ybr_transcripts, temp2, by = colnames(temp2)[-8], all.x = T)
    ybr_transcripts$orf[ybr_transcripts$name %in% temp$name] <-
      paste(ybr_transcripts$orf[ybr_transcripts$name %in% temp$name], annot$GeneID[i], sep = ',')
  }
}

save.image(file = sprintf('%s/tx_archi.RData',out_path))


##### READ ARCHI FILE
load(file = sprintf('%s/tx_archi.RData',out_path))
ybr_transcripts$orf <- str_remove(ybr_transcripts$orf, 'NA,')
head(ybr_transcripts)

temp <- ybr_transcripts %>%
  group_by(sample, orf) %>%
  count() %>% data.frame()
tx_counts <- data.frame(tx = unique(ybr_transcripts$orf))
col.names <- NULL
for (s in unique(temp$sample)) {
  tx_counts <- merge(tx_counts, temp[temp$sample == s, c('orf','n')],
                     by.x = 'tx', by.y = 'orf', all.x = T,
                     suffixes = c('',sprintf('_%s',s)))
  col.names <- c(col.names, s)
}
colnames(tx_counts) <- c('tx',col.names)
tx_counts$n_orfs <- str_count(tx_counts$tx, ',') + 1
head(tx_counts)


tx_counts %>%
  ggplot(aes(x = n_orfs)) +
  geom_line(stat = 'density') +
  scale_x_log10(breaks = c(seq(1,10,1),100))

samples %>%
  ggplot(aes(x = stop - start)) +
  geom_line(stat = 'density') +
  facet_wrap(.~sample) +
  scale_x_log10()


tx_counts %>%
  ggplot(aes(x = syn2, y = stop2)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

##### YBR and Friends
hello <- ybr_transcripts[(str_detect(ybr_transcripts$orf, 'orf14869,') | ybr_transcript$orf == 'orf14869') & !is.na(ybr_transcripts$orf),]

hi <- tx_counts[str_detect(tx_counts$tx, 'orf14869,') | tx_counts$tx == 'orf14869',]
head(hi)
sum(rowSums(hi[,c(2:5)], na.rm = T))
