##### LOAD ALL OVERLAP RESULTS
overlap_gene.one_two_all <- read.csv(file = 'output/pull_down_one_and_two_all_samples_overlap.csv')
overlap_gene.one_two_mng <- read.csv(file = 'output/pull_down_one_and_two_mNG_samples_overlap.csv')
overlap_gene.two_ha_mng <- read.csv(file = 'output/pull_down_two_HA_and_mNG_samples_overlap.csv')
overlap_gene.two_ha <- read.csv(file = 'output/pull_down_two_HA_samples_overlap.csv')
overlap_gene.two_mng <- read.csv(file = 'output/pull_down_two_mNG_samples_overlap.csv')

gcn1_pi <- read.csv(file = 'input/ms_data/GCN1_physical_interactions.csv')
gcn4_targets <- read.csv(file = 'input/GCN4_targets.csv', col.names = c('DBID','orf_name','organism','Gene','description'))
gcn4_degs <- read.csv(file = 'input/GCN4_DEGs.csv')

###
colnames(overlap_gene.two_mng)
overlap_gene.two_mng %>%
  filter(test_overlap == 1, peptide_cnts_control == 0) %>%
  nrow()

overlap_gene.two_mng %>%
  filter(norm_lfc > 0)

hi <- overlap_gene.two_ha_mng %>%
  filter(lfc > 1, test_overlap == 1
         # ,
         # Gene %notin% unique(gcn1_pi$Interactor.1)
         )
write.csv(hi, file = 'output/pull_down_enrichment_two_ha_mng.csv', row.names = F)

overlap_gene.two_ha_mng %>%
  filter(Gene == 'GDT1')


### DEA AND COEXP GCN1
## MON2, ARL1
coexp.ybr %>%
  filter(GENENAME %in% c('MON2','ARL1'))

dea_results %>%
  filter(orf_name %in% c('YNL297C','YBR164C'), !is.na(DE),
         contrast == 'YBR3_BY4741_DEL_vs_YBR3_BY4741_WT')

## GCN pathway
coexp.ybr %>%
  filter(GENENAME %in% c('GCN1','GCN2','GCN3','GCN4','GCN20','SPF1'))

dea_results %>%
  filter(orf_name %in% c('YGL195W','YDR283C','YKR026C','YEL009C','YFR009W','YEL031W','YBR196C','YCR059C'), !is.na(DE),
         # abs(lfc) > log2(1),
         contrast == 'YBR3_BY4741_DEL_vs_YBR3_BY4741_WT') #%>%
  #group_by(orf_name, DE) %>%
  #count()

## IRE1 HAC1
dea_results %>%
  filter(orf_name %in% c('YHR079C','YFL031W','YJR132W'),
         # abs(lfc) > log2(1),
         contrast == 'YBR3_BY4741_DEL_vs_YBR3_BY4741_WT')

dea_results %>%
  filter(orf_name %in% unique(gcn4_targets$orf_name), !is.na(DE), abs(lfc) > 1,
         contrast == 'YBR3_BY4741_DEL_vs_YBR3_BY4741_WT') %>% nrow()
dea_results %>%
  filter(!is.na(DE), abs(lfc) > 1,
         contrast == 'YBR3_BY4741_DEL_vs_YBR3_BY4741_WT') %>% nrow()


##### ANNOTATED COEXP
coexp.ybr.annotates <- coexp.ybr %>% filter(!is.na(orf_name))
row.names(coexp.ybr.annotates) <- NULL

write.csv(coexp.ybr.annotates, file = 'output/ydl_coexpression_network2.csv', row.names = F)

##### HA SAMPLES
pddat %>%
  filter(ybr_tag %in% c('HA','YBR_HA'), protein_present == 1) %>%
  group_by(Attempt, background, ybr_tag, replicate, Gene) %>%
  count() %>%
  group_by(Attempt, background, ybr_tag, Gene) %>%
  count() %>% filter(n == 3) %>%
  group_by(Gene) %>%
  count() %>% filter(n == 2)





