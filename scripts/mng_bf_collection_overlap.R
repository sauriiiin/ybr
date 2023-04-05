

mng.col <- read.csv(file = 'input/mNG_collection_init.csv')
mng.col <- mng.col[!str_detect(mng.col$Sequence.Name, '_RV'),]
mng.col <- data.frame(orf_name = unique(mng.col))
head(mng.col)

bfcol <- dbGetQuery(conn, 'select * from BARFLEX_SPACE_AGAR_180313
                    union
                    select * from PROTOGENE_COLLECTION')
head(bfcol)


mng.col <- merge(mng.col, bfcol, by = 'orf_name')
write.csv(mng.col, file = 'input/mNG_BF_overlap_for_nelsio.csv', row.names = FALSE)
