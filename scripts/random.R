
hello <- read.csv(file = 'input/zapregulon.csv', header = T)
head(hello)
dbWriteTable(conn, 'WANG2018_ZAPREGULON', hello[,c(-64,-65)])


hello <- read.csv(file = 'input/znproteome.csv', header = T)
head(hello)
dbWriteTable(conn, 'WANG2018_ZINCPROTEOME', hello)


hello <- read.csv(file = 'input/co_wt_vs_atg_degs.csv', header = T)
head(hello)
dbWriteTable(conn, 'RNASEQ_YBRATG_COMMON_DEA', hello)


hello <- read.table(file = 'input/ZAP1_targets_high_throughput.txt', header = T, sep = '\t', skip = 8)
head(hello)
dbWriteTable(conn, 'SGD_HT_ZAPREGULON', hello[,c(1:5,14,15)])


hello <- read.table(file = 'https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2011-12-12-r125/MediaObjects/13059_2011_2751_MOESM2_ESM.TXT', sep = '', fill = T)
dbWriteTable(conn, 'SCER_TF', 
             data.frame(TF = str_to_upper(unique(str_split(hello$V1[str_detect(hello$V1, '_')], '_', simplify = T)[,1]))))


hello <- read.csv(file = 'input/scertfs.csv')
dbWriteTable(conn, 'SCER_TF', hello, overwrite = T)
