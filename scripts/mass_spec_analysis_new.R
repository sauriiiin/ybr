

# min pep 1	min prob 90	min crap expt 15

library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(UpSetR)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(RMariaDB)

`%notin%` <- Negate(`%in%`)

out_path <- "~/R/Projects/ybr/output"
fig_path <- "~/R/Projects/ybr/output/figs"

source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

GFP_crampone <- read_excel("input/interactions/crapome.xlsx")

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

##### LOADING DATA
sample_info <- read.csv(file = 'input/ms_data/221225_sample_info.csv')
# sample_info <- sample_info[,-7]
sample_info$Sample.ID.number <- str_replace_all(sample_info$Sample.ID.number, '-', '.')
sample_info$Sample.description <- str_remove_all(sample_info$Sample.description, ' ')
head(sample_info)

hatag_yeastonly <- read.csv(file = 'input/ms_data/221225_YBRPullDown_HATag_YeastonlyDB.csv')
hatag_yeastplus <- read.csv(file = 'input/ms_data/221225_YBRPullDown_HATag_YeastDB_OpenProt.csv')
mngtag_yeastonly <- read.csv(file = 'input/ms_data/221225_YBRPullDown_mNGTag_YeastonlyDB.csv')
mngtag_yeastplus <- read.csv(file = 'input/ms_data/221225_YBRPullDown_mNGTag_YeastDB_OpenProt.csv')

head(mngtag_yeastplus)
colnames(mngtag_yeastplus)

pulldown2 <- rbind(data.frame(tag = 'mNG', protein_db = 'yeast', merge(sample_info, melt(mngtag_yeastonly, id.vars = colnames(mngtag_yeastonly)[c(1:12,45)],
                                                                            variable.name = 'Sample.ID.number', value.name = 'spectral_counts'), by = 'Sample.ID.number')),
      data.frame(tag = 'mNG', protein_db = 'yeast_openprot', merge(sample_info, melt(mngtag_yeastplus, id.vars = colnames(mngtag_yeastplus)[c(1:12,45)],
                                                                                     variable.name = 'Sample.ID.number', value.name = 'spectral_counts'), by = 'Sample.ID.number')),
      data.frame(tag = 'HA', protein_db = 'yeast', merge(sample_info, melt(hatag_yeastonly, id.vars = colnames(hatag_yeastonly)[c(1:12,61)],
                                                                           variable.name = 'Sample.ID.number', value.name = 'spectral_counts'), by = 'Sample.ID.number')),
      data.frame(tag = 'HA', protein_db = 'yeast_openprot', merge(sample_info, melt(hatag_yeastplus, id.vars = colnames(hatag_yeastplus)[c(1:12,61)],
                                                                                    variable.name = 'Sample.ID.number', value.name = 'spectral_counts'), by = 'Sample.ID.number')))
colnames(pulldown2)

pulldown2 %>%
  ggplot(aes(x = Protein.Probability)) +
  geom_density() +
  facet_wrap(tag ~ protein_db)

GFP_crampone %>%
  ggplot(aes(x = NUM_EXPT)) +
  geom_density()

pulldown2 %>%
  filter(Protein.Probability > 0.9) %>%
  # filter(Gene == 'YBR196C-A') %>%
  # filter(Gene %notin% unique(GFP_crampone$GENE[GFP_crampone$NUM_EXPT >= 15])) %>%
  group_by(tag, protein_db, Sample.description, Protein, Protein.ID, Entry.Name, Gene, Protein.Length, Organism) %>%
  summarise(spectral_counts = sum(spectral_counts, na.rm = T), .groups = 'keep') %>%
  filter(spectral_counts > 0) %>%
  group_by(tag, protein_db, Sample.description) %>%
  count() %>% data.frame()

pulldown2 %>%
  filter(Protein.Probability > 0.9) %>%
  # filter(Gene == 'YBR196C-A') %>%
  # filter(Gene %notin% unique(GFP_crampone$GENE[GFP_crampone$NUM_EXPT >= 15])) %>%
  group_by(tag, protein_db, Sample.description, Protein, Protein.ID, Entry.Name, Gene, Protein.Length, Organism) %>%
  summarise(spectral_counts = sum(spectral_counts, na.rm = T), .groups = 'keep') %>%
  filter(spectral_counts > 100) %>%
  data.frame() 

pulldown2 %>%
  filter(Protein.Probability > 0.9) %>%
  # filter(Gene == 'YBR196C-A') %>%
  # filter(Gene %notin% unique(GFP_crampone$GENE[GFP_crampone$NUM_EXPT >= 15])) %>%
  group_by(tag, protein_db, Sample.description, Protein, Protein.ID, Entry.Name, Gene, Protein.Length, Organism) %>%
  summarise(spectral_counts = sum(spectral_counts, na.rm = T), .groups = 'keep') %>%
  filter(spectral_counts > 100) %>%
  ggplot(aes(x = spectral_counts)) +
  geom_histogram() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(tag ~ protein_db * Sample.description)

temp <- pulldown2 %>%
  # mutate(Sample.description = case_when(str_detect(Sample.description, 'YBR-HA') ~ 'HA_yeast_YBR-HA',
  #                                       str_detect(Sample.description, 'untagged') ~ 'HA_yeast_untagged',
  #                                       TRUE ~ Sample.description)) %>%
  filter(protein_db == 'yeast', Protein.Probability > 0.9) %>%
  group_by(tag, protein_db, Sample.description, Protein, Protein.ID, Entry.Name, Gene, Protein.Length, Organism) %>%
  summarise(spectral_counts = sum(spectral_counts, na.rm = T), .groups = 'keep')


temp <- temp %>%
  mutate(protein_presence = case_when(spectral_counts >  0 ~ 1,
                                      spectral_counts == 0 ~ 0))
head(temp)

temp$Gene[str_detect(temp$Gene, ',')]

pulldown3 <- data.frame(Gene = unique(temp$Gene))
colsss <- NULL
for (t in unique(temp$tag)) {
  for (db in unique(temp$protein_db[temp$tag == t])) {
    for (s in unique(temp$Sample.description[temp$tag == t & temp$protein_db == db])) {
      colsss <- c(colsss, paste(t,db,s,sep = '_'))
      
      pulldown3 <- merge(pulldown3,
            temp[temp$tag == t & temp$protein_db == db & temp$Sample.description == s, c('Gene', 'protein_presence')],
            by = 'Gene', all.x = T)
      
      
    }
  }
}
pulldown3 <- pulldown3[c(-1,-2),]
colnames(pulldown3) <- c('Gene',colsss)
pulldown3[is.na(pulldown3)] <- 0
head(pulldown3)


upset(pulldown3, nsets = 10, order.by = 'freq')#, nintersects = 18)

upset(pulldown3[,c('Gene', colnames(pulldown3)[str_detect(colnames(pulldown3), 'mNG')])], nsets = 4,
      order.by = 'freq', nintersects = 18)


hello <- pulldown3[,c('Gene', colnames(pulldown3)[str_detect(colnames(pulldown3), 'mNG')])] %>%
  mutate(control = `mNG_yeast_mNG-doa10C` + `mNG_yeast_mNG-WTC`,
         ybr = `mNG_yeast_YBR-mNG-doa10C` + `mNG_yeast_YBRmNG-WTC`) %>%
  filter(control == 0, ybr > 0)

hi <- bitr(hello$Gene, fromType = "GENENAME",
           toType = c("ENTREZID","GENENAME","ENSEMBL"),
           OrgDb = org.Sc.sgd.db)
paste(hi$ENSEMBL, collapse = ', ')





upset(pulldown3[,c('Gene', colnames(pulldown3)[str_detect(colnames(pulldown3), 'HA')])], nsets = 6,
      order.by = 'freq')

upset(pulldown3[,c('Gene', colnames(pulldown3)[str_detect(colnames(pulldown3), 'HA')])] %>%
  mutate(untagged = HA_yeast_untaggedyeastrep1 + HA_yeast_untaggedyeastrep2 + HA_yeast_untaggedyeastrep3) %>%
  filter(untagged == 0), nsets = 6, order.by = 'freq')

hello <- pulldown3[,c('Gene', colnames(pulldown3)[str_detect(colnames(pulldown3), 'HA')])] %>%
  mutate(untagged = HA_yeast_untaggedyeastrep1 + HA_yeast_untaggedyeastrep2 + HA_yeast_untaggedyeastrep3,
         tagged = `HA_yeast_YBR-HAyeastrep1` + `HA_yeast_YBR-HAyeastrep2` + `HA_yeast_YBR-HAyeastrep3`) %>%
  filter(untagged == 0, tagged > 0)

hi <- bitr(hello$Gene, fromType = "GENENAME",
     toType = c("ENTREZID","GENENAME","ENSEMBL"),
     OrgDb = org.Sc.sgd.db)
paste(hi$ENSEMBL, collapse = ', ')
