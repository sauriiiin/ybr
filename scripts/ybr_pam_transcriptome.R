##### YBR196C-A SHORT READ (ATTEMPT2) TRANSCRIPTOME
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 11/15/2022

source(file = 'scripts/initialize.R')

##### COMPILE COUNT DATA
load('~/R/Projects/ybr/output/YBR2_featureCount_raw_sgd_withoverlap.RData')
raw_cnts <- as.data.frame(fc$counts)
colnames(raw_cnts) <- paste('YBR2', c(paste('WT',seq(1,3),sep='_'),
                                      paste('PAM1',seq(1,3),sep='_'),
                                      paste('PAM2',seq(1,3),sep='_'),
                                      paste('ATG',seq(1,3),sep='_'),
                                      paste('DEL',seq(1,3),sep='_'),
                                      paste('OE',seq(1,3),sep='_')), sep = '_')
raw_cnts <- raw_cnts[str_detect(rownames(raw_cnts), 'Y') &
                       !str_detect(rownames(raw_cnts), 'tY'),]
head(raw_cnts)


##### PCA
pca_obj = prcomp(raw_cnts)
pca <- as.data.frame(pca_obj[2]$rotation)
pca[,"samples"] <- str_split(row.names(pca), '_', simplify = T)[,2]
pca_sum <- summary(pca_obj)

pca %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(aes(color=samples), size = 3) +
  labs(title = 'Uncorrected Counts',
       x = sprintf('PC1 (%0.2f%%)', pca_sum$importance[2,1] * 100),
       y = sprintf('PC2 (%0.2f%%)', pca_sum$importance[2,2] * 100)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom')


##### MEAN COUNTS
raw_cnts_sum <- raw_cnts %>%
  mutate(orf_name = row.names(.)) %>%
  melt(id.vars = 'orf_name', variable.name = 'sid', value.name = 'count')
temp <- cbind(unique(as.character(raw_cnts_sum$sid)), str_split(unique(raw_cnts_sum$sid), '_', simplify = T)) %>% data.frame()
colnames(temp) <- c('sid','batch','sample','replicate')
raw_cnts_sum <- merge(temp, raw_cnts_sum, by = 'sid', all = T)
raw_cnts_sum <- raw_cnts_sum %>%
  group_by(batch, sample, orf_name) %>%
  summarize(mean_cnt = mean(count, na.rm = T),
            median_cnt = median(count, na.rm = T),
            std_cnt = sd(count, na.rm = T),
            .groups = 'keep') %>%
  mutate(cv_cnt = std_cnt/mean_cnt) %>%
  data.frame()
head(raw_cnts_sum)

sample_cmbn <- t(combn(unique(raw_cnts_sum$sample), 2))

raw_cnts_sum_mat <- NULL
for (i in seq(1,nrow(sample_cmbn))) {
  temp <- merge(raw_cnts_sum[,c('sample','orf_name','median_cnt')] %>% filter(sample == sample_cmbn[i,1]),
        raw_cnts_sum[,c('sample','orf_name','median_cnt')] %>% filter(sample == sample_cmbn[i,2]),
        by = 'orf_name', suffixes = c('_1','_2'))
  raw_cnts_sum_mat <- rbind(raw_cnts_sum_mat, temp)
}
raw_cnts_sum_mat <- data.frame(raw_cnts_sum_mat)

raw_cnts_sum_mat %>%
  ggplot(aes(x = median_cnt_1, y = median_cnt_2)) +
  geom_point() +
  stat_cor(method = 'spearman') +
  facet_wrap(sample_1 ~ sample_2)


##### DEA
ybr2_wt_samples = paste(paste('YBR2','WT',sep='_'), seq(1,3), sep = '_')
ybr2_pam1_samples = paste(paste('YBR2','PAM1',sep='_'), seq(1,3), sep = '_')
ybr2_pam2_samples = paste(paste('YBR2','PAM2',sep='_'), seq(1,3), sep = '_')
ybr2_atg_samples = paste(paste('YBR2','ATG',sep='_'), seq(1,3), sep = '_')
ybr2_del_samples = paste(paste('YBR2','DEL',sep='_'), seq(1,3), sep = '_')
ybr2_oe_samples = paste(paste('YBR2','OE',sep='_'), seq(1,3), sep = '_')

ybr_proplus_samples = c(ybr2_wt_samples, ybr2_pam1_samples, ybr2_pam2_samples)
ybr_prominus_samples = c(ybr2_atg_samples, ybr2_del_samples)

ybr2_dea_results <- NULL
for (i in seq(1,nrow(sample_cmbn))) {
  temp <- run_edgeR(data=raw_cnts,
                    group_a_name=sample_cmbn[i,2], group_a_samples=paste(paste('YBR2',sample_cmbn[i,2],sep='_'), seq(1,3), sep = '_'),
                    group_b_name=sample_cmbn[i,1], group_b_samples=paste(paste('YBR2',sample_cmbn[i,1],sep='_'), seq(1,3), sep = '_'))
  ybr2_dea_results <- rbind(ybr2_dea_results, temp)
}
temp <- run_edgeR(data=raw_cnts,
                  group_a_name='PRO_PLUS', group_a_samples=ybr_proplus_samples,
                  group_b_name='PRO_MINUS', group_b_samples=ybr_prominus_samples)
ybr2_dea_results <- rbind(ybr2_dea_results, temp)
temp <- run_edgeR(data=raw_cnts,
                  group_a_name='PAM', group_a_samples=c(ybr2_pam1_samples,ybr2_pam2_samples),
                  group_b_name='ATG', group_b_samples=ybr2_atg_samples)
ybr2_dea_results <- rbind(ybr2_dea_results, temp)
temp <- run_edgeR(data=raw_cnts,
                  group_a_name='PAM', group_a_samples=c(ybr2_pam1_samples,ybr2_pam2_samples),
                  group_b_name='DEL', group_b_samples=ybr2_del_samples)
ybr2_dea_results <- rbind(ybr2_dea_results, temp)
temp <- run_edgeR(data=raw_cnts,
                  group_a_name='PAM', group_a_samples=c(ybr2_pam1_samples,ybr2_pam2_samples),
                  group_b_name='ATG_DEL', group_b_samples=c(ybr2_atg_samples, ybr2_del_samples))
ybr2_dea_results <- rbind(ybr2_dea_results, temp)
temp <- run_edgeR(data=raw_cnts,
                  group_a_name='WT', group_a_samples=ybr2_wt_samples,
                  group_b_name='ATG_DEL', group_b_samples=c(ybr2_atg_samples, ybr2_del_samples))
ybr2_dea_results <- rbind(ybr2_dea_results, temp)

ybr2_dea_results$DE[ybr2_dea_results$fdr <= 0.05 & ybr2_dea_results$lfc > 0] <- 'Up'
ybr2_dea_results$DE[ybr2_dea_results$fdr <= 0.05 & ybr2_dea_results$lfc < 0] <- 'Down'

ybr2_dea_results %>%
  filter(!is.na(DE)) %>%
  group_by(contrast, DE) %>%
  count() %>% data.frame()

dbWriteTable(conn, 'RNASEQ_YBRPAM_DEA', ybr2_dea_results, overwrite = T)


