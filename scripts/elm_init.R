

arm_protein_doms <- read.csv(file = 'input/ms_data/HA_36/arm_proteins_all_domains.csv',
                             col.names = c('orf_name','std_name','org_name','domain_name','domain_description','start','end'))

elm_int <- readr::read_tsv(file = 'input/elm_analysis/elm_interactions.tsv') %>% data.frame()
head(elm_int)

elm_int <- elm_int %>%
  dplyr::filter(taxonomyElm %in% c('559292(Saccharomyces cerevisiae S288c)', '4932(Saccharomyces cerevisiae)'),
                taxonomyDomain %in% c('559292(Saccharomyces cerevisiae S288c)', '4932(Saccharomyces cerevisiae)'))

elm_int %>%
  dplyr::filter(Elm %in% c('CLV_PCSK_KEX2_1',
                           'DOC_CYCLIN_yCln2_LP_2',
                           'DOC_MAPK_MEF2A_6',
                           'LIG_Pex14_2',
                           'LIG_SUMO_SIM_par_1',
                           'TRG_ER_diArg_1',
                           'DOC_MAPK_DCC_7',
                           'LIG_WD40_WDR5_VDV_2',
                           'MOD_Plk_1',
                           'TRG_ENDOCYTIC_2')) %>%
  group_by(Elm, Domain) %>% count()

elm_int %>%
  dplyr::filter(Elm %in% c('CLV_PCSK_KEX2_1',
                           'DOC_CYCLIN_yCln2_LP_2',
                           'DOC_MAPK_MEF2A_6',
                           'LIG_Pex14_2',
                           'LIG_SUMO_SIM_par_1',
                           'TRG_ER_diArg_1',
                           'DOC_MAPK_DCC_7',
                           'LIG_WD40_WDR5_VDV_2',
                           'MOD_Plk_1',
                           'TRG_ENDOCYTIC_2'),
                Domain %in% unique(arm_protein_doms$domain_name))



elm_classes <- readr::read_tsv(file = 'input/elm_analysis/elm_classes.tsv') %>% data.frame()
ybr_friends_elm <- read.csv(file = 'input/elm_analysis/ybr_friends_elm.csv')

merge(ybr_friends_elm, elm_int, by.x = 'elm_name', by.y = 'Elm') %>%
  group_by(orf_name, elm_name, Domain) %>% count()

write.csv(merge(ybr_friends_elm %>%
  group_by(orf_name, elm_name) %>%
  count() %>%
  group_by(elm_name) %>% count() %>%
  dplyr::filter(n > 0), elm_classes[,c(2:4)], by.x = 'elm_name', by.y = 'ELMIdentifier'),
  file = 'output/ybr_friends_elm_classes.csv',
  row.names = F)

