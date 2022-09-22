

#import sample
sample_30a<-read.delim('/home/aar75/rna_seq/longreads/nanopore_ybr/20220819/30a.bed',header=F,col.names=c('chr','start','stop','name','score','strand'))
sample_30b<-read.delim('/home/aar75/rna_seq/longreads/nanopore_ybr/20220819/30b.bed',header=F,col.names=c('chr','start','stop','name','score','strand'))
sample_30c<-read.delim('/home/aar75/rna_seq/longreads/nanopore_ybr/20220819/30c.bed',header=F,col.names=c('chr','start','stop','name','score','strand'))
sample_30d<-read.delim('/home/aar75/rna_seq/longreads/nanopore_ybr/20220819/30d.bed',header=F,col.names=c('chr','start','stop','name','score','strand'))

samples <- rbind(cbind(sample_30b, sample = 'wt'),
                 cbind(sample_30a, sample = 'crispy'),
                 cbind(sample_30d, sample = 'syn'),
                 cbind(sample_30c, sample = 'stop')) %>% data.frame()
head(samples)
samples %>%
  mutate(n_bases = stop - start + 1) %>%
  group_by(sample) %>%
  # summarize(basecalled = sum(n_bases))
  count()


#import annotation info for all orfs you are interested in quantifying 
orf_info<-readRDS('/home/aar75/rna_seq/Salmon_4_6_22/translatome.rds')
head(orf_info)

#loop through all orfs in annotation 
annot$GeneID[i]

annot[i,]
i <- 14086 #ybr
counts <- data.frame(annot, wt = 0, crispy = 0, syn = 0, stop = 0)
# samples_matrix <- samples
for(i in 1:nrow(counts)){
  temp <- samples[(samples$chr == counts$Chr[i] & samples$strand == counts$Strand[i]) &
               ((samples$start <= counts$Start[i] & samples$stop >= counts$Start[i]) |
                  (samples$start <= counts$End[i] & samples$stop >= counts$End[i]) |
                  (samples$start >= counts$Start[i] & samples$stop <= counts$End[i]) |
                  (samples$start <= counts$Start[i] & samples$stop >= counts$End[i])),]
  if (nrow(temp) > 0) {
    # if (i %in% seq(14086-10,14086+10,1)) {
    #   temp2 <- cbind(temp, hello = 1)
    #   colnames(temp2) <- c(colnames(temp2)[-8], counts$GeneID[i])
    #   samples_matrix <- merge(samples_matrix, temp2, by = colnames(temp2)[-8], all.x = T)
    # }
    # 
    temp <- temp %>%
      group_by(sample) %>% 
      count() %>%
      data.frame()
    if (length(temp$n[temp$sample == 'wt']) > 0) {counts$wt[i] <- temp$n[temp$sample == 'wt']}
    if (length(temp$n[temp$sample == 'crispy']) > 0) {counts$crispy[i] <- temp$n[temp$sample == 'crispy']}
    if (length(temp$n[temp$sample == 'syn']) > 0) {counts$syn[i] <- temp$n[temp$sample == 'syn']}
    if (length(temp$n[temp$sample == 'stop']) > 0) {counts$stop[i] <- temp$n[temp$sample == 'stop']}
  }
}
head(counts)
# head(samples_matrix)

reads <- merge(orfid_crossmap, counts, by.x = 'gene_id', by.y = 'GeneID', all.y = T)
reads <- merge(orfs, reads, by = 'orf_id', all.y = T)
head(reads)


##### TRANSCRIPT ARCHITECHTURE FOR YBR LOCUS
loci_range <- seq(14086-1,14086+4,1)
ybr_transcripts <- samples
for(i in loci_range){
  temp <- samples[(samples$chr == annot$Chr[i] & samples$strand == annot$Strand[i]) &
                    ((samples$start <= annot$Start[i] & samples$stop >= annot$Start[i]) |
                       (samples$start <= annot$End[i] & samples$stop >= annot$End[i]) |
                       (samples$start >= annot$Start[i] & samples$stop <= annot$End[i]) |
                       (samples$start <= annot$Start[i] & samples$stop >= annot$End[i])),]
  if (nrow(temp) > 0) {
    temp2 <- cbind(temp, hello = 1)
    colnames(temp2) <- c(colnames(temp2)[-8], annot$GeneID[i])
    ybr_transcripts <- merge(ybr_transcripts, temp2, by = colnames(temp2)[-8], all.x = T)
  }
}
head(ybr_transcripts)
ybr_transcripts <- ybr_transcripts[rowSums(ybr_transcripts[,8:ncol(ybr_transcripts)], na.rm = T) > 0,]
ybr_transcripts$t_length <- ybr_transcripts$stop - ybr_transcripts$start
head(ybr_transcripts)


# order the transcripts based on start site and length
ybr_transcripts <- ybr_transcripts[order(ybr_transcripts$sample, ybr_transcripts$stop, -ybr_transcripts$t_length),]
ybr_transcripts <- ybr_transcripts %>%
  mutate(t_order = seq(1,nrow(ybr_transcripts),1))

ybr_transcripts %>%
  group_by(sample) %>% count()

ybr_transcripts[,c('start','stop','name','strand','sample')] %>%
  melt(id.vars = c('strand','sample','name')) %>%
  summarize(quantiles_0_50_100 = quantile(value, c(0,0.5,1)))
  
ybr_loci <- annot[loci_range,]
ybr_loci <- merge(merge(ybr_loci, orfid_crossmap, by.x = 'GeneID', by.y = 'gene_id'),
                  orf_info[,c('orf_id','transcript','gene_systematic_name','rfc','is_transient')],
                  by = 'orf_id')
ybr_loci <- ybr_loci[order(ybr_loci$Start),]
ybr_loci <- ybr_loci %>%
  mutate(o_order = seq(1,nrow(ybr_loci),1),
         o_label_position = (Start+End)/2,
         o_label = sprintf('%s (%s)', GeneID, gene_systematic_name))

loci_zoom <- c(612000, 616000)

plot.transcript.crispy <- ybr_transcripts[,c('start','stop','name','strand','sample','t_length','t_order')] %>%
  melt(id.vars = c('strand','sample','name','t_order','t_length')) %>%
  filter(sample == 'stop') %>%
  ggplot(aes(x = value, y = as.factor(t_order))) +
  geom_vline(xintercept = unique(c(ybr_loci$Start, ybr_loci$End)), linetype = 'dashed', col = 'blue', size = 0.3) +
  geom_line() +
  facet_grid(.~sample) +
  coord_cartesian(xlim = loci_zoom) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

plot.ybr.loci <- ybr_loci %>%
  melt(id.vars = colnames(ybr_loci)[c(-4,-5)]) %>%
  ggplot(aes(x = value, y = o_order)) +
  geom_vline(xintercept = unique(c(ybr_loci$Start, ybr_loci$End)), linetype = 'dashed', col = 'blue', size = 0.3) +
  geom_line(aes(group = o_order)) +
  geom_text(aes(x = o_label_position, y = o_order+0.3, label = o_label), size = 2, hjust = 0.5) +
  # scale_y_continuous(breaks = seq(1,20,1)) +
  coord_cartesian(xlim = loci_zoom) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())


plot_grid(plot.transcript.crispy, plot.ybr.loci,
          nrow = 2, rel_heights = c(5,1),
          align = 'v')
  

