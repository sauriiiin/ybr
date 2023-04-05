
library(readxl)
library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)

path <- "/home/sbp29/R/Projects/ybr/input/ms_data/221225_YBR_MS_LaneWise_TagSplit.xlsx"
sheets <- excel_sheets(path = path)
sheets <- sheets[!str_detect(sheets, 'Open')]

##### CHARACTERIZE MASS SPEC DATA

hatag_pep <- read_xlsx(path = path, sheet = sheets[1]) %>% data.frame()
hatag_pro <- read_xlsx(path = path, sheet = sheets[2]) %>% data.frame()
mngtag_pep <- read_xlsx(path = path, sheet = sheets[3]) %>% data.frame()
mngtag_pro <- read_xlsx(path = path, sheet = sheets[4]) %>% data.frame()


hatag_bp <- barplot(colSums(hatag_pep[9:14]),
              xlab = "Sample", 
              ylab = "Total Peptide Counts",
              main = "HA Tag")
text(hatag_bp, 1000, colSums(hatag_pep[9:14]))


mngtag_bp <- barplot(colSums(mngtag_pep[9:12]),
                    xlab = "Sample", 
                    ylab = "Total Peptide Counts",
                    main = "mNG Tag")
text(mngtag_bp, 2000, colSums(mngtag_pep[9:12]))
colSums(mngtag_pep[9:12])


##### COMMON
hatag_pro[,c('Gene','Protein.Probability','WT_rep1', 'WT_rep2', 'WT_rep3', 'YBR_HA_rep1', 'YBR_HA_rep2', 'YBR_HA_rep3')] %>%
  filter(`Protein.Probability` > 0.9) %>%
  melt(id.vars = c('Gene','Protein.Probability')) %>%
  filter(value > 0, variable %in% c('YBR_HA_rep1', 'YBR_HA_rep2')) %>%
  group_by(Gene) %>%
  count() %>%
  filter(n == 2) %>% nrow()


mngtag_pro[,c('Gene','Protein.Probability','WT_mNG', 'doa10_mNG', 'WT_YBR_mNG', 'doa10_YBR_mNG')] %>%
  filter(`Protein.Probability` > 0.9) %>%
  melt(id.vars = c('Gene','Protein.Probability')) %>%
  filter(value > 0, variable %in% c('WT_YBR_mNG', 'doa10_YBR_mNG')) %>%
  group_by(Gene) %>%
  count() %>%
  filter(n == 2) %>% nrow()



