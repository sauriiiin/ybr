

coldata <- read.csv('~/RNASeq/YBR/fastq_files/YBR1/samples.csv', stringsAsFactors = FALSE)
colnames(coldata) <- c('filename','name')
# coldata$name <- factor(coldata$name, levels = c("WT_plas","CRISPY","WT","DEL"))
coldata$name <- factor(coldata$name, levels = c("WT","KAN","CRISPY","CRISPY-ATG"))
coldata$replicate <- as.factor(rep(c(1,2,3),4))
coldata$sample <- paste(coldata$name, coldata$replicate, sep = '_')
rownames(coldata) <- coldata$sample

coldata1 <- coldata
coldata2 <- coldata

head(coldata1)
head(coldata2)


load('~/R/Projects/ybr/output/YBR1_featureCount_raw.RData')
cnts1 <- as.data.frame(fc$counts)
colnames(cnts1) <- paste('YBR1', coldata1$sample, sep = '_')
head(cnts1)

load('~/R/Projects/ybr/output/YBR3_featureCount_raw.RData')
cnts2 <- as.data.frame(fc$counts)
colnames(cnts2) <- paste('YBR3', coldata2$sample, sep = '_')
head(cnts2)


coldata <- data.frame(rbind(cbind(coldata1, attempt = 'YBR1'),
                            cbind(coldata2, attempt = 'YBR3')))
coldata$sample_name <- paste(coldata$attempt, coldata$sample, sep = '_')

sum(rownames(cnts1) != rownames(cnts2))

cnts <- cbind(cnts1, cnts2)

##### FILTERING
cds <- read.fasta(file = "~/RNASeq/YBR/scer_cds.fa")
cds <- str_remove(names(cds), '_mRNA')
cnts <- cnts[rownames(cnts) %in% cds,]
sprintf('There are %d genes with more than 10 reads in total.' ,sum(rowSums(cnts) >= 10))
sprintf('There are %d genes with more than 10 reads/sample.' ,sum(rowSums(cnts < 10) == 0))

keep <- rowSums(cnts) >= 10
# keep <- rowSums(cnts < 10) == 0
cnts <- cnts[keep,]

##### RAW CORRELATION
raw.cor <-cor(cnts)
data.frame(raw.cor)
col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue","#00007F"))
corrplot(raw.cor,
         type="lower", order="original",
         col=col1(10),
         method = 'circle',
         is.corr = FALSE,
         cl.lim = c(round(min(raw.cor),1),1))

##### NORMALIZING COUNTS
rownames(coldata) <- NULL
dds <- DESeqDataSetFromMatrix(countData = as.matrix(cnts),
                              colData = coldata %>% mutate(sampss = paste(attempt, name, sep = '_')),
                              design = ~sampss)
dds <- DESeq(dds)

##### SAMPLE DISTANCE
vsd <- vst(dds)
sampleDists <- dist(t(assay(vsd)))
# cor(assay(vsd), method = 'pearson')

sampleDistMatrix <- as.matrix( sampleDists )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         show_rownames = T,
         show_colnames = T)

#### PCA
plotPCA(vsd, intgroup = c("name")) +
  geom_text(aes(label = name), position = position_nudge(y = 2))



##### CORRELATION PLOTS
cnts_nrm <- counts(dds, normalized=TRUE)
head(cnts_nrm)


cnts %>%
  data.frame() %>%
  ggplot(aes(x = YBR1_WT_1, y = YBR3_WT_plas_1)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  coord_cartesian(xlim = c(0,200000),
                  ylim = c(0,200000))

##### DEA
res <- data.frame()
for (c in resultsNames(dds)[2:length(resultsNames(dds))]){
  # cat(sprintf('Log fold change, DE, GO and KEGG results for %s.\n',str_remove(c, 'name_')))
  cat(sprintf('Processing %s...\n',str_remove(c, 'name_')))
  temp <- lfcShrink(dds, coef=c, type="apeglm")
  DESeq2::plotMA(temp, ylim = c(-5, 5), main=str_remove(c, 'name_'))
  cat(sprintf('When comparing %s there are %d genes with more than 1 log fold change at FDR of 5%%. With %d upregulated and %d downregulated.\n\n',
              str_remove(c, 'name_'),
              sum(temp$padj < 0.05 & abs(temp$log2FoldChange) >= 1, na.rm = T),
              sum(temp$padj < 0.05 & temp$log2FoldChange >= 1, na.rm = T),
              sum(temp$padj < 0.05 & temp$log2FoldChange <= -1, na.rm = T)))
  temp$orf_name <- rownames(temp)
  rownames(temp) <- NULL
  res <- rbind(res, data.frame(temp, contrast = str_remove(c, 'name_')))
}


###### RES1 VS RES3
res1 <- read.csv(file = '~/R/Projects/ybr/output/ybr1_res.csv')
res2 <- read.csv(file = '~/R/Projects/ybr/output/ybr3_res.csv')

unique(res1$contrast)
unique(res2$contrast)

merge(res1 %>% filter(contrast == "CRISPY.ATG_vs_WT"), 
      res2 %>% filter(contrast == "CRISPY_vs_WT_plas"),
      by = 'orf_name', all = T, suffixes = c('_ybr1','_ybr3')) %>%
  ggplot(aes(x = log2FoldChange_ybr1, y = log2FoldChange_ybr3)) +
  geom_point() +
  stat_cor()


