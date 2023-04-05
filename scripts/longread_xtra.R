library(dplyr)

sequencing_summary <- rbind(
  cbind(read.table(file = '/home/sbp29/RNASeq/YBR/minION/2208/sample_30b/basecalling/sequencing_summary.txt', header = T), sample = 'wt2'),
  cbind(read.table(file = '/home/sbp29/RNASeq/YBR/minION/2208/sample_30a_/basecalling/sequencing_summary.txt', header = T), sample = 'crispy2'),
  cbind(read.table(file = '/home/sbp29/RNASeq/YBR/minION/2208/sample_30d/basecalling/sequencing_summary.txt', header = T), sample = 'syn2'),
  cbind(read.table(file = '/home/sbp29/RNASeq/YBR/minION/2208/sample_30c/basecalling/sequencing_summary.txt', header = T), sample = 'stop2')
) %>%
  data.frame()
head(sequencing_summary)

sequencing_report <- merge(merge(merge(sequencing_summary %>%
              group_by(sample, passes_filtering) %>% count() %>% data.frame(),
            sequencing_summary %>%
              group_by(sample) %>% count() %>% data.frame(),
            by = 'sample', suffixes = c('','_total')) %>%
        mutate(filter_perc = n/n_total * 100),
      sequencing_summary %>%
        group_by(sample, passes_filtering) %>%
        summarize(min_qscore = min(mean_qscore_template),
          median_qscore = median(mean_qscore_template),
                  max_qscore = max(mean_qscore_template),
                  .groups = 'keep'),
      by = c('sample','passes_filtering')),
      samples %>%
        group_by(sample) %>%
        count() %>% data.frame(), by = 'sample', suffixes = c('','_aligned')) %>%
  mutate(n_not_aligned = n - n_aligned,
         align_pass_perc = n_aligned/n * 100,
         align_perc = n_aligned/n_total * 100)
write.csv(sequencing_report, file = 'output/minion_report.csv')

samples2 <- merge(samples, sequencing_summary, by.x = 'name', by.y = 'read_id', all.x = T)

length(unique(sequencing_summary$read_id[sequencing_summary$sample == 'crispy2']))

head(sample_30a2)

merge(sample_30a2, sequencing_summary[,c('filename','read_id','passes_filtering','mean_qscore_template')], by.x = 'name', by.y = 'read_id') %>%
  # group_by(passes_filtering) %>% count()
  # summarize(q = mean(mean_qscore_template, na.rm = T))
  ggplot(aes(x = mean_qscore_template)) +
  geom_line(stat = 'density', trim = T)






loci_zoom <- c(612000, 616000)


plot.loci <- ybr_loci %>%
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

# plot.transcripts.wt <- ybr_transcripts[,c('start','stop','name','strand','sample','t_length','t_order')] %>%
#   melt(id.vars = c('strand','sample','name','t_order','t_length')) %>%
#   filter(sample == 'wt') %>%
#   ggplot(aes(x = value, y = as.factor(t_order))) +
#   geom_vline(xintercept = unique(c(ybr_loci$Start, ybr_loci$End)), linetype = 'dashed', col = 'blue', size = 0.3) +
#   geom_line(size = 0.3) +
#   facet_grid(.~sample) +
#   coord_cartesian(xlim = loci_zoom) +
#   theme_minimal() +
#   theme(panel.grid = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank())
# plot.ybr.loci.wt <- plot_grid(plot.transcripts.wt, plot.loci,
#                            nrow = 2, rel_heights = c(5,1),
#                            align = 'v')


plot.transcripts.wt2 <- ybr_transcripts[,c('start','stop','name','strand','sample','t_length','t_order')] %>%
  melt(id.vars = c('strand','sample','name','t_order','t_length')) %>%
  filter(sample == 'wt2') %>%
  ggplot(aes(x = value, y = as.factor(t_order))) +
  geom_vline(xintercept = unique(c(ybr_loci$Start, ybr_loci$End)), linetype = 'dashed', col = 'blue', size = 0.3) +
  geom_line(size = 0.3) +
  facet_grid(.~sample) +
  coord_cartesian(xlim = loci_zoom) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
plot.ybr.loci.wt2 <- plot_grid(plot.transcripts.wt2, plot.loci,
                              nrow = 2, rel_heights = c(5,1),
                              align = 'v')


# plot.transcripts.crispy <- ybr_transcripts[,c('start','stop','name','strand','sample','t_length','t_order')] %>%
#   melt(id.vars = c('strand','sample','name','t_order','t_length')) %>%
#   filter(sample == 'crispy') %>%
#   ggplot(aes(x = value, y = as.factor(t_order))) +
#   geom_vline(xintercept = unique(c(ybr_loci$Start, ybr_loci$End)), linetype = 'dashed', col = 'blue', size = 0.3) +
#   geom_line(size = 0.3) +
#   facet_grid(.~sample) +
#   coord_cartesian(xlim = loci_zoom) +
#   theme_minimal() +
#   theme(panel.grid = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank())
# plot.ybr.loci.crispy <- plot_grid(plot.transcripts.crispy, plot.loci,
#                               nrow = 2, rel_heights = c(5,1),
#                               align = 'v')

plot.transcripts.crispy2 <- ybr_transcripts[,c('start','stop','name','strand','sample','t_length','t_order')] %>%
  melt(id.vars = c('strand','sample','name','t_order','t_length')) %>%
  filter(sample == 'crispy2') %>%
  ggplot(aes(x = value, y = as.factor(t_order))) +
  geom_vline(xintercept = unique(c(ybr_loci$Start, ybr_loci$End)), linetype = 'dashed', col = 'blue', size = 0.3) +
  geom_line(size = 0.3) +
  facet_grid(.~sample) +
  coord_cartesian(xlim = loci_zoom) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
plot.ybr.loci.crispy2 <- plot_grid(plot.transcripts.crispy2, plot.loci,
                                  nrow = 2, rel_heights = c(5,1),
                                  align = 'v')

# plot.transcripts.syn <- ybr_transcripts[,c('start','stop','name','strand','sample','t_length','t_order')] %>%
#   melt(id.vars = c('strand','sample','name','t_order','t_length')) %>%
#   filter(sample == 'syn') %>%
#   ggplot(aes(x = value, y = as.factor(t_order))) +
#   geom_vline(xintercept = unique(c(ybr_loci$Start, ybr_loci$End)), linetype = 'dashed', col = 'blue', size = 0.3) +
#   geom_line(size = 0.3) +
#   facet_grid(.~sample) +
#   coord_cartesian(xlim = loci_zoom) +
#   theme_minimal() +
#   theme(panel.grid = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank())
# plot.ybr.loci.syn <- plot_grid(plot.transcripts.syn, plot.loci,
#                                    nrow = 2, rel_heights = c(5,1),
#                                    align = 'v')

plot.transcripts.syn2 <- ybr_transcripts[,c('start','stop','name','strand','sample','t_length','t_order')] %>%
  melt(id.vars = c('strand','sample','name','t_order','t_length')) %>%
  filter(sample == 'syn2') %>%
  ggplot(aes(x = value, y = as.factor(t_order))) +
  geom_vline(xintercept = unique(c(ybr_loci$Start, ybr_loci$End)), linetype = 'dashed', col = 'blue', size = 0.3) +
  geom_line(size = 0.3) +
  facet_grid(.~sample) +
  coord_cartesian(xlim = loci_zoom) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
plot.ybr.loci.syn2 <- plot_grid(plot.transcripts.syn2, plot.loci,
                               nrow = 2, rel_heights = c(5,1),
                               align = 'v')

# plot.transcripts.stop <- ybr_transcripts[,c('start','stop','name','strand','sample','t_length','t_order')] %>%
#   melt(id.vars = c('strand','sample','name','t_order','t_length')) %>%
#   filter(sample == 'stop') %>%
#   ggplot(aes(x = value, y = as.factor(t_order))) +
#   geom_vline(xintercept = unique(c(ybr_loci$Start, ybr_loci$End)), linetype = 'dashed', col = 'blue', size = 0.3) +
#   geom_line(size = 0.3) +
#   facet_grid(.~sample) +
#   coord_cartesian(xlim = loci_zoom) +
#   theme_minimal() +
#   theme(panel.grid = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank())
# plot.ybr.loci.stop <- plot_grid(plot.transcripts.stop, plot.loci,
#                                nrow = 2, rel_heights = c(5,1),
#                                align = 'v')


plot.transcripts.stop2 <- ybr_transcripts[,c('start','stop','name','strand','sample','t_length','t_order')] %>%
  melt(id.vars = c('strand','sample','name','t_order','t_length')) %>%
  filter(sample == 'stop2') %>%
  ggplot(aes(x = value, y = as.factor(t_order))) +
  geom_vline(xintercept = unique(c(ybr_loci$Start, ybr_loci$End)), linetype = 'dashed', col = 'blue', size = 0.3) +
  geom_line(size = 0.3) +
  facet_grid(.~sample) +
  coord_cartesian(xlim = loci_zoom) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
plot.ybr.loci.stop2 <- plot_grid(plot.transcripts.stop2, plot.loci,
                               nrow = 2, rel_heights = c(5,1),
                               align = 'v')



# plot.ybr.loci <- plot_grid(plot.ybr.loci.wt, plot.ybr.loci.wt2,
#           plot.ybr.loci.crispy, plot.ybr.loci.crispy2,
#           plot.ybr.loci.syn, plot.ybr.loci.syn2,
#           plot.ybr.loci.stop, plot.ybr.loci.stop2,
#           nrow = 4, ncol = 2)
# 
# ggsave(sprintf("%s/YBR_LOCUS_TRANSCRIPTS.jpg",fig_path), plot.ybr.loci,
#        height = 1000, width = 500, units = 'mm',
#        bg = 'white',
#        dpi = 300,
#        limitsize = FALSE)


plot.ybr.loci <- plot_grid(plot.ybr.loci.wt2,
                           plot.ybr.loci.crispy2,
                           plot.ybr.loci.syn2,
                           plot.ybr.loci.stop2,
                           nrow = 2, ncol = 2)

ggsave(sprintf("%s/YBR_LOCUS_TRANSCRIPTS2.jpg",fig_path), plot.ybr.loci,
       height = 500, width = 400, units = 'mm',
       bg = 'white',
       dpi = 300,
       limitsize = FALSE)


###### MATRIX PLOT
ybr_transcripts %>%
  filter(orf14869 == 1)


###### COUNTING
# i <- 14086
loci_range <- seq(1,nrow(annot),1) #seq(14086-1,14086+4,1)
min_overlap <- 10
ybr_transcripts <- samples
ybr_transcripts$orfs <- NA
for(i in loci_range){
  temp <- samples[(samples$chr == annot$Chr[i] & samples$strand == annot$Strand[i]) &
                    (((samples$start <= annot$Start[i] & samples$stop >= annot$Start[i]) & (samples$stop - annot$Start[i]) > min_overlap) |
                       ((samples$start <= annot$End[i] & samples$stop >= annot$End[i]) & (annot$End[i] - samples$start) > min_overlap) |
                       ((samples$start >= annot$Start[i] & samples$stop <= annot$End[i]) & (samples$start - samples$stop) > min_overlap) |
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
# ybr_transcripts <- ybr_transcripts[rowSums(ybr_transcripts[,8:ncol(ybr_transcripts)], na.rm = T) > 0,]
ybr_transcripts$t_length <- ybr_transcripts$stop - ybr_transcripts$start

ybr_transcripts <- ybr_transcripts[order(ybr_transcripts$sample, ybr_transcripts$stop, -ybr_transcripts$t_length),]
ybr_transcripts <- ybr_transcripts %>%
  mutate(t_order = seq(1,nrow(ybr_transcripts),1))
head(ybr_transcripts)

ybr_transcripts %>%
  group_by(sample) %>% count()

ybr_loci <- annot[loci_range,]
ybr_loci <- merge(merge(ybr_loci, orfid_crossmap, by.x = 'GeneID', by.y = 'gene_id'),
                  orf_info[,c('orf_id','transcript','gene_systematic_name','rfc','is_transient')],
                  by = 'orf_id')
ybr_loci <- ybr_loci[order(ybr_loci$Start),]
ybr_loci <- ybr_loci %>%
  mutate(o_order = seq(1,nrow(ybr_loci),1),
         o_label_position = (Start+End)/2,
         o_label = sprintf('%s (%s)', GeneID, gene_systematic_name))





# ?matrix
# 
# hello <- matrix(data = 0, nrow = nrow(ybr_transcripts), ncol = nrow(annot),
#                 dimnames = list(annot$GeneID, ybr_transcripts$sample_id))
