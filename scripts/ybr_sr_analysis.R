##### YBR196C-A ALL SHORT READ TRANSCRIPTOME
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 11/15/2022

source(file = 'scripts/initialize.R')

##### GATHER DATA
sgd_zap1_regulon <- read.table('input/ZAP1_targets.txt', sep = '\t', skip = 8, header = T)
sgd_zap1_regulon <- sgd_zap1_regulon[,c(1:5,14,15)]
head(sgd_zap1_regulon)
length(unique(sgd_zap1_regulon$Target.Systematic.Name))

sgd_zap1_ht_regulon <- read.table('input/ZAP1_targets_high_throughput.txt', sep = '\t', skip = 8, header = T)
sgd_zap1_ht_regulon <- sgd_zap1_ht_regulon[,c(1:5,14,15)]
head(sgd_zap1_ht_regulon)
length(unique(sgd_zap1_ht_regulon$Target.Systematic.Name))

zn_proteome <- read.csv('input/znproteome.csv')
head(zn_proteome)
length(unique(zn_proteome$Systematic.Name))

all_genes <- dbGetQuery(conn, 'select * from RNASEQ_YBRATG_ALL_GENES')
all_coip <- dbGetQuery(conn, 'select * from YBR_COIP_RESULTS')


##### COMPILE COUNT DATA
load('~/R/Projects/ybr/output/YBR1_featureCount_raw_sgd.RData')
cnts1 <- as.data.frame(fc$counts[,c(1:3,7:12)])
colnames(cnts1) <- paste('YBR1', paste('BY4742', c(paste('WT',seq(1,3),sep='_'),paste('DEL',seq(1,3),sep='_'),paste('ATG',seq(1,3),sep='_')), sep = '_'), sep = '_')
head(cnts1)

load('~/R/Projects/ybr/output/YBR2_featureCount_raw_sgd.RData')
cnts2 <- as.data.frame(fc$counts)
colnames(cnts2) <- paste('YBR2', paste('BY4742', c(paste('WT',seq(1,3),sep='_'),
                                      paste('PAM1',seq(1,3),sep='_'),
                                      paste('PAM2',seq(1,3),sep='_'),
                                      paste('ATG',seq(1,3),sep='_'),
                                      paste('DEL',seq(1,3),sep='_'),
                                      paste('OE',seq(1,3),sep='_')), sep = '_'),
                         sep = '_')
head(cnts2)

load('~/R/Projects/ybr/output/YBR3_featureCount_raw_sgd.RData')
cnts3 <- as.data.frame(fc$counts)
colnames(cnts3) <- paste('YBR3', c(paste('BY4742', c(paste('WT',seq(1,3),sep='_'),paste('ATG',seq(1,3),sep='_')), sep = '_'),
                                   paste('BY4741', c(paste('WT',seq(1,3),sep='_'),paste('DEL',seq(1,3),sep='_')), sep = '_')),
                         sep = '_')
head(cnts3)

raw_cnts <- cbind(cnts1, cnts2, cnts3)
raw_cnts <- raw_cnts[str_detect(rownames(raw_cnts), 'Y') &
                       !str_detect(rownames(raw_cnts), 'tY'),]

raw_cnts[row.names(raw_cnts) %in% c('YLR303W','YEL021W'),]

##### BATCH CORRECTION
batch_group = as.numeric(as.factor(str_split(colnames(raw_cnts), '_', simplify = T)[,1]))
bg_group = as.numeric(as.factor(str_split(colnames(raw_cnts), '_', simplify = T)[,2]))
crisp_group = c(rep(1,9),rep(2,18),rep(1,6),rep(3,6))
sample_group = as.numeric(as.factor(str_split(colnames(raw_cnts), '_', simplify = T)[,3]))
covariate_matrix = cbind(bg_group, sample_group)

pca_uncorrected_obj = prcomp(raw_cnts[,c(1:9,28:33)])
pca_uncorrected <- as.data.frame(pca_uncorrected_obj[2]$rotation)
pca_uncorrected[,"batch"] <- str_split(colnames(raw_cnts[,c(1:9,28:33)]), '_', simplify = T)[,1]
pca_uncorrected[,"bg"] <- str_split(colnames(raw_cnts[,c(1:9,28:33)]), '_', simplify = T)[,2]
pca_uncorrected[,"samples"] <- str_split(colnames(raw_cnts[,c(1:9,28:33)]), '_', simplify = T)[,3]
pca_uncorrected_sum <- summary(pca_uncorrected_obj)

# corrected_cnts = ComBat_seq(counts = as.matrix(raw_cnts), batch = batch_group, group = NULL, covar_mod = covariate_matrix)
# corrected_cnts = ComBat_seq(counts = as.matrix(raw_cnts), batch = batch_group, group = bg_group, full_mod = T)
corrected_cnts = ComBat_seq(counts = as.matrix(raw_cnts[,c(1:9,28:33)]), batch = batch_group[c(1:9,28:33)], group = NULL, full_mod = T) #sample_group[c(1:9,28:33)]
pca_corrected_obj = prcomp(corrected_cnts)
pca_corrected <- as.data.frame(pca_corrected_obj[2]$rotation)
pca_corrected[,"batch"] <- str_split(colnames(raw_cnts[,c(1:9,28:33)]), '_', simplify = T)[,1]
pca_corrected[,"bg"] <- str_split(colnames(raw_cnts[,c(1:9,28:33)]), '_', simplify = T)[,2]
pca_corrected[,"samples"] <- str_split(colnames(raw_cnts[,c(1:9,28:33)]), '_', simplify = T)[,3]
pca_corrected_sum <- summary(pca_corrected_obj)

plot.uncorrected.pca <- pca_uncorrected %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(aes(color=batch, shape=samples), size = 3) +
  labs(title = 'Uncorrected Counts',
       x = sprintf('PC1 (%0.2f%%)', pca_uncorrected_sum$importance[2,1] * 100),
       y = sprintf('PC2 (%0.2f%%)', pca_uncorrected_sum$importance[2,2] * 100)) +
  facet_wrap(.~bg) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom')

plot.corrected.pca <- pca_corrected %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(aes(color=batch, shape=samples), size = 3) +
  labs(title = 'Corrected Counts',
       x = sprintf('PC1 (%0.2f%%)', pca_corrected_sum$importance[2,1] * 100),
       y = sprintf('PC2 (%0.2f%%)', pca_corrected_sum$importance[2,2] * 100)) +
  facet_wrap(.~bg) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom')

plot.samples.pca <- ggarrange(plot.uncorrected.pca, plot.corrected.pca,
                              nrow = 2, align = 'hv',
                              common.legend = T, legend = 'bottom')
plot.samples.pca

##### MERGE THE RAW AND CORRECTED CNTS
corrected_cnts <- cbind(corrected_cnts, raw_cnts[,c(10:27,34:39)])

##### DEA
sample_cmbn <- t(combn(unique(str_remove(str_remove(str_remove(colnames(corrected_cnts), '_1'), '_2'), '_3')), 2))
sample_cmbn <- sample_cmbn[((str_detect(sample_cmbn[,1], 'YBR1') & str_detect(sample_cmbn[,2], 'YBR1')) |
               (str_detect(sample_cmbn[,1], 'YBR2') & str_detect(sample_cmbn[,2], 'YBR2')) |
               (str_detect(sample_cmbn[,1], 'YBR3') & str_detect(sample_cmbn[,2], 'YBR3'))) &
              ((str_detect(sample_cmbn[,1], 'BY4741') & str_detect(sample_cmbn[,2], 'BY4741')) |
                 (str_detect(sample_cmbn[,1], 'BY4742') & str_detect(sample_cmbn[,2], 'BY4742'))),]

dea_results <- NULL
for (i in seq(1,nrow(sample_cmbn))) {
  temp <- run_edgeR(data=corrected_cnts,
                    group_a_name=sample_cmbn[i,1], group_a_samples=paste(sample_cmbn[i,1], seq(1,3), sep = '_'),
                    group_b_name=sample_cmbn[i,2], group_b_samples=paste(sample_cmbn[i,2], seq(1,3), sep = '_'))
  temp$batch = str_split(sample_cmbn[i,1], '_', simplify = T)[,1]
  temp$background = str_split(sample_cmbn[i,1], '_', simplify = T)[,2]
  temp$comparison = sprintf('%s_vs_%s',str_split(sample_cmbn[i,2], '_', simplify = T)[,3], str_split(sample_cmbn[i,1], '_', simplify = T)[,3])
  dea_results <- rbind(dea_results, temp)
}
head(dea_results)

temp <- run_edgeR(data=corrected_cnts,
                  group_a_name='YBR1_YBR3_BY4742_WT', group_a_samples=c(paste('YBR1_BY4742_WT', seq(1,3), sep = '_'),
                                                                        paste('YBR3_BY4742_WT', seq(1,3), sep = '_')),
                  group_b_name='YBR1_YBR3_BY4742_ATG', group_b_samples=c(paste('YBR1_BY4742_ATG', seq(1,3), sep = '_'),
                                                                         paste('YBR3_BY4742_ATG', seq(1,3), sep = '_')))
temp$batch = 'YBR1/3'
temp$background = 'BY4742'
temp$comparison = 'ATG_vs_WT'
dea_results <- rbind(dea_results, temp)

temp <- run_edgeR(data=corrected_cnts,
                  group_a_name='YBR1_YBR2_YBR3_BY4742_WT', group_a_samples=c(paste('YBR1_BY4742_WT', seq(1,3), sep = '_'),
                                                                             paste('YBR2_BY4742_WT', seq(1,3), sep = '_'),
                                                                             paste('YBR3_BY4742_WT', seq(1,3), sep = '_')),
                  group_b_name='YBR1_YBR2_YBR3_BY4742_ATG', group_b_samples=c(paste('YBR1_BY4742_ATG', seq(1,3), sep = '_'),
                                                                              paste('YBR2_BY4742_ATG', seq(1,3), sep = '_'),
                                                                              paste('YBR3_BY4742_ATG', seq(1,3), sep = '_')))
temp$batch = 'YBR1/2/3'
temp$background = 'BY4742'
temp$comparison = 'ATG_vs_WT'
dea_results <- rbind(dea_results, temp)

# temp <- run_edgeR(data=corrected_cnts,
#                   group_a_name='YBR1_YBR3_BY4742_BY4741_WT', group_a_samples=c(paste('YBR1_BY4742_WT', seq(1,3), sep = '_'),
#                                                                         paste('YBR3_BY4741_WT', seq(1,3), sep = '_')),
#                   group_b_name='YBR1_YBR3_BY4742_BY4741_DEL', group_b_samples=c(paste('YBR1_BY4742_DEL', seq(1,3), sep = '_'),
#                                                                          paste('YBR3_BY4741_DEL', seq(1,3), sep = '_')))
# temp$batch = 'YBR1/3'
# temp$background = 'BY4741/2'
# temp$comparison = 'DEL_vs_WT'
# dea_results <- rbind(dea_results, temp)
# 
# temp <- run_edgeR(data=corrected_cnts,
#                   group_a_name='YBR2_YBR3_BY4741_WT', group_a_samples=c(paste('YBR2_BY4741_WT', seq(1,3), sep = '_'),
#                                                                                paste('YBR3_BY4741_WT', seq(1,3), sep = '_')),
#                   group_b_name='YBR2_YBR3_BY4741_DEL', group_b_samples=c(paste('YBR2_BY4741_DEL', seq(1,3), sep = '_'),
#                                                                                 paste('YBR3_BY4741_DEL', seq(1,3), sep = '_')))
# temp$batch = 'YBR2/3'
# temp$background = 'BY4741'
# temp$comparison = 'DEL_vs_WT'
# dea_results <- rbind(dea_results, temp)

dea_results$DE[dea_results$fdr <= 0.05 & dea_results$lfc > 0] <- 'Up'
dea_results$DE[dea_results$fdr <= 0.05 & dea_results$lfc < 0] <- 'Down'

dea_results %>%
  filter(orf_name == 'YJL056C', comparison %in% c('DEL_vs_WT', 'ATG_vs_WT'))

dbWriteTable(conn, 'RNASEQ_YBR_DEA', dea_results, overwrite = T)

##### DE COUNTS
dea_cnts <- merge(dea_results %>%
  filter(!is.na(DE), comparison %in% c('DEL_vs_WT', 'ATG_vs_WT')) %>%
  group_by(batch, background, comparison) %>%
  count(),
  dea_results %>%
    filter(!is.na(DE), comparison %in% c('DEL_vs_WT', 'ATG_vs_WT')) %>%
    group_by(batch, background, comparison, DE) %>%
    count(), by = c('batch','background','comparison'), all = T, suffixes = c('_total_DE','_DE'))

##### ZINC REGULON, ZAP PROTEINS & COIP COUNTS
dea_cnts <- merge(dea_cnts,merge(merge(dea_results %>%
                                         filter(orf_name %in% unique(sgd_zap1_ht_regulon$Target.Systematic.Name), comparison %in% c('DEL_vs_WT', 'ATG_vs_WT')) %>%
                                         group_by(batch, background, comparison) %>%
                                         count(),
                                       dea_results %>%
                                         filter(orf_name %in% unique(sgd_zap1_ht_regulon$Target.Systematic.Name), !is.na(DE), comparison %in% c('DEL_vs_WT', 'ATG_vs_WT')) %>%
                                         group_by(batch, background, comparison) %>%
                                         count(), by = c('batch','background','comparison'), all = T, suffixes = c('_total_zap1_regulon_ht', '_total_DE_zap1_regulon_ht')),
                                 dea_results %>%
                                   filter(orf_name %in% unique(sgd_zap1_ht_regulon$Target.Systematic.Name), !is.na(DE), comparison %in% c('DEL_vs_WT', 'ATG_vs_WT')) %>%
                                   group_by(batch, background, comparison,DE) %>%
                                   count(), by = c('batch','background','comparison'), all = T, suffixes = c('', '_DE_zap1_regulon_ht')),
                  by = c('batch','background','comparison','DE'), all = T)

dea_cnts <- merge(dea_cnts,merge(merge(dea_results %>%
              filter(orf_name %in% unique(sgd_zap1_regulon$Target.Systematic.Name), comparison %in% c('DEL_vs_WT', 'ATG_vs_WT')) %>%
              group_by(batch, background, comparison) %>%
              count(),
            dea_results %>%
              filter(orf_name %in% unique(sgd_zap1_regulon$Target.Systematic.Name), !is.na(DE), comparison %in% c('DEL_vs_WT', 'ATG_vs_WT')) %>%
              group_by(batch, background, comparison) %>%
              count(), by = c('batch','background','comparison'), all = T, suffixes = c('_total_zap1_regulon', '_total_DE_zap1_regulon')),
      dea_results %>%
        filter(orf_name %in% unique(sgd_zap1_regulon$Target.Systematic.Name), !is.na(DE), comparison %in% c('DEL_vs_WT', 'ATG_vs_WT')) %>%
        group_by(batch, background, comparison,DE) %>%
        count(), by = c('batch','background','comparison'), all = T, suffixes = c('', '_DE_zap1_regulon')),
      by = c('batch','background','comparison','DE'), all = T, suffixes = c('_DE_zap1_regulon_ht','_DE_zap1_regulon'))



dea_cnts <- merge(dea_cnts,merge(merge(dea_results %>%
                             filter(orf_name %in% unique(zn_proteome$Systematic.Name), comparison %in% c('DEL_vs_WT', 'ATG_vs_WT')) %>%
                             group_by(batch, background, comparison) %>%
                             count(),
                           dea_results %>%
                             filter(orf_name %in% unique(zn_proteome$Systematic.Name), !is.na(DE), comparison %in% c('DEL_vs_WT', 'ATG_vs_WT')) %>%
                             group_by(batch, background, comparison) %>%
                             count(), by = c('batch','background','comparison'), all = T, suffixes = c('_total_zn_proteome', '_total_DE_zn_proteome')),
                     dea_results %>%
                       filter(orf_name %in% unique(zn_proteome$Systematic.Name), !is.na(DE), comparison %in% c('DEL_vs_WT', 'ATG_vs_WT')) %>%
                       group_by(batch, background, comparison,DE) %>%
                       count(), by = c('batch','background','comparison'), all = T, suffixes = c('', '_DE_zinc_proteome')),
      by = c('batch','background','comparison','DE'), all = T)


dea_cnts <- merge(dea_cnts,merge(merge(dea_results %>%
                             filter(orf_name %in% unique(all_coip$ORF), comparison %in% c('DEL_vs_WT', 'ATG_vs_WT')) %>%
                             group_by(batch, background, comparison) %>%
                             count(),
                           dea_results %>%
                             filter(orf_name %in% unique(all_coip$ORF), !is.na(DE), comparison %in% c('DEL_vs_WT', 'ATG_vs_WT')) %>%
                             group_by(batch, background, comparison) %>%
                             count(), by = c('batch','background','comparison'), all = T, suffixes = c('_total_ybr_coip', '_total_DE_ybr_coip')),
                     dea_results %>%
                       filter(orf_name %in% unique(all_coip$ORF), !is.na(DE), comparison %in% c('DEL_vs_WT', 'ATG_vs_WT')) %>%
                       group_by(batch, background, comparison,DE) %>%
                       count(), by = c('batch','background','comparison'), all = T, suffixes = c('', '_DE_ybr_coip')),
      by = c('batch','background','comparison','DE'), all = T, suffixes = c('_DE_zinc_proteome', '_DE_ybr_coip'))

dea_cnts$n_total_genes <- length(unique(dea_results$orf_name))

##### ZINC REGULON, ZAP PROTEINS & COIP ENRICHMENT
head(dea_cnts)

for (batch in unique(dea_cnts$batch)) {
  for (bg in unique(dea_cnts$background[dea_cnts$batch == batch])) {
    for (comp in unique(dea_cnts$comparison[dea_cnts$batch == batch & dea_cnts$background == bg])) {
      temp <- dea_cnts[dea_cnts$batch == batch & dea_cnts$background == bg & dea_cnts$comparison == comp,]
      
      a <- mean(temp$n_total_DE_zap1_regulon_ht, na.rm = T)
      b <- mean(temp$n_total_DE, na.rm = T) - a
      c <- mean(temp$n_total_zap1_regulon_ht, na.rm = T) - a
      d <- mean(temp$n_total_genes, na.rm = T) - (a + b + c)
      temp_mat <- matrix(c(a, b, c, d),
                         nrow = 2,
                         dimnames = list(ZAP1=c("yes","no"),
                                         DE=c("yes","no")))
      temp_en <- fisher.test(temp_mat, alternative = "greater")
      dea_cnts$p_total_DE_zap1_regulon_ht[dea_cnts$batch == batch & dea_cnts$background == bg & dea_cnts$comparison == comp] <- temp_en$p.value
      dea_cnts$or_total_DE_zap1_regulon_ht[dea_cnts$batch == batch & dea_cnts$background == bg & dea_cnts$comparison == comp] <- temp_en$estimate[[1]]
      
      
      a <- mean(temp$n_total_DE_zap1_regulon, na.rm = T)
      b <- mean(temp$n_total_DE, na.rm = T) - a
      c <- mean(temp$n_total_zap1_regulon, na.rm = T) - a
      d <- mean(temp$n_total_genes, na.rm = T) - (a + b + c)
      temp_mat <- matrix(c(a, b, c, d),
                         nrow = 2,
                         dimnames = list(ZAP1=c("yes","no"),
                                         DE=c("yes","no")))
      temp_en <- fisher.test(temp_mat, alternative = "greater")
      dea_cnts$p_total_DE_zap1_regulon[dea_cnts$batch == batch & dea_cnts$background == bg & dea_cnts$comparison == comp] <- temp_en$p.value
      dea_cnts$or_total_DE_zap1_regulon[dea_cnts$batch == batch & dea_cnts$background == bg & dea_cnts$comparison == comp] <- temp_en$estimate[[1]]
      
      
      a <- mean(temp$n_total_DE_zn_proteome, na.rm = T)
      b <- mean(temp$n_total_DE, na.rm = T) - a
      c <- mean(temp$n_total_zn_proteome, na.rm = T) - a
      d <- mean(temp$n_total_genes, na.rm = T) - (a + b + c)
      temp_mat <- matrix(c(a, b, c, d),
                         nrow = 2,
                         dimnames = list(ZN=c("yes","no"),
                                         DE=c("yes","no")))
      temp_en <- fisher.test(temp_mat, alternative = "greater")
      dea_cnts$p_total_DE_zn_proteome[dea_cnts$batch == batch & dea_cnts$background == bg & dea_cnts$comparison == comp] <- temp_en$p.value
      dea_cnts$or_total_DE_zn_proteome[dea_cnts$batch == batch & dea_cnts$background == bg & dea_cnts$comparison == comp] <- temp_en$estimate[[1]]
    }
  }
}

dbWriteTable(conn, 'RNASEQ_YBR_DEA_COUNTS', dea_cnts, overwrite = T)

head(zn_proteome)

##### ATPASE
head(all_genes)
merge(dea_results, all_genes, by.x = 'orf_name', by.y = 'ORF') %>%
  filter(!is.na(DE), str_detect(DESCRIPTION, 'V-ATPase'), comparison %in% c('DEL_vs_WT', 'ATG_vs_WT')) %>%
  filter(orf_name %in% all_coip$ORF) %>%
  group_by(batch, background, comparison) %>%
  # group_by(GENENAME) %>%
  count() %>% data.frame()

sum(str_detect(all_genes$DESCRIPTION, 'V-ATPase'))


##### IRON COPPER
head(all_genes)
merge(dea_results, all_genes, by.x = 'orf_name', by.y = 'ORF') %>%
  filter(!is.na(DE), str_detect(DESCRIPTION, 'cu'), comparison %in% c('DEL_vs_WT', 'ATG_vs_WT')) %>%
  group_by(batch, background, comparison) %>%
  count() %>% data.frame()

sum(str_detect(all_genes$DESCRIPTION, 'cu'))

