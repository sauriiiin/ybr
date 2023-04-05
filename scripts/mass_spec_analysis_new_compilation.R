
library(readxl)

sample_info <- read.csv(file = '/home/sbp29/R/Projects/ybr/input/ms_data/221225_sample_info.csv')
sample_info$Sample.ID.number <- str_replace_all(sample_info$Sample.ID.number, '-', '_')
sample_info <- sample_info[1:80,]
sample_info$Sample.ID.number[sample_info$Sample.ID.number == 'AFO_9'] <- 'AFO_09'

path <- "input/ms_data/122522 O'Donnell Tag Comparison_YeastDB + OpenProt.xlsx"
sheets <- excel_sheets(path = path)

# str_detect(colnames(dat.in), 'AFO_9\\b')


for (sheet.name in sheets) {
  dat.in <- read_xlsx(path = path, sheet = sheet.name, skip = 6) %>% data.frame()
  colnames(dat.in)[str_detect(colnames(dat.in), 'AFO_')] <-
    c(str_extract_all(colnames(dat.in)[str_detect(colnames(dat.in), 'AFO_')], 'AFO_\\d+', simplify = T))
  if (sum(str_detect(colnames(dat.in), 'AFO_9\\b')) > 0) {
    colnames(dat.in)[str_detect(colnames(dat.in), 'AFO_9\\b')] <- 'AFO_09'
  }
  
  dat.out <- dat.in[,!str_detect(colnames(dat.in), 'AFO_')]
  for (s in unique(sample_info$Sample.description)) {
    for (r in unique(sample_info$replicate[sample_info$Sample.description == s])) {
      ids <- sample_info$Sample.ID.number[sample_info$Sample.description == s & sample_info$replicate == r]
      if (sum(ids %in% colnames(dat.in)) == length(ids)) {
        dat.out[sprintf('%s',s)] <- rowSums(dat.in[,ids])
        write.csv(dat.out, file = sprintf('input/ms_data/122522_YBR_MS_YeastOpenProtDB_%s.csv',sheet.name), row.names = F)
      }
    }
  }
}

##### COMBINING HA AND MNG TAG DATA
path <- "input/ms_data/221225_YBR_MS_LaneWise_TagSplit.xlsx"
sheets <- excel_sheets(path = path)[5:8]


dat.hatag.pep <- read_xlsx(path = path, sheet = sheets[1]) %>% data.frame()
dat.hatag.pro <- read_xlsx(path = path, sheet = sheets[2]) %>% data.frame()
dat.mngtag.pep <- read_xlsx(path = path, sheet = sheets[3]) %>% data.frame()
dat.mngtag.pro <- read_xlsx(path = path, sheet = sheets[4]) %>% data.frame()

colnames(dat.hatag.pep) <- str_remove_all(colnames(dat.hatag.pep), "\\.")
colnames(dat.hatag.pro) <- str_remove_all(colnames(dat.hatag.pro), "\\.")
colnames(dat.mngtag.pep) <- str_remove_all(colnames(dat.mngtag.pep), "\\.")
colnames(dat.mngtag.pro) <- str_remove_all(colnames(dat.mngtag.pro), "\\.")

dat.pep <- merge(dat.hatag.pep[,-4], dat.mngtag.pep[,-4], by = c("Sequence", "ChargeStates", "Gene", "Protein", "ProteinID", "ProteinDescription"),
      all = T, suffixes = c('_HA_Tag','_mNG_Tag'))

colnames(dat.mngtag.pro)
colnames(dat.hatag.pro)

dat.pro <- merge(dat.hatag.pro, dat.mngtag.pro, by = c("Protein", "ProteinID", "EntryName", "Gene", "ProteinLength", "Organism"),
      all = T, suffixes = c('_HA_Tag','_mNG_Tag'))










