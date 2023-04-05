`%notin%` <- Negate(`%in%`)

annot <- read.table(file = '/home/sbp29/RNASeq/YBR/translated_orfs_with_transient_status.gff3')
annot <- annot[,c(9,1,4,5,7)]
annot$V9 <- str_remove(str_split(annot$V9, ';', simplify = T)[,1], 'ID=')
colnames(annot) <- c('GeneID','Chr','Start','End','Strand')
annot$Chr <- str_remove(annot$Chr, 'chr')


fc1 <- featureCounts(files="~/RNASeq/YBR/fastq_files/YBR3/rawdata/quant_star/3SP1_quantAligned.sortedByCoord.out.bam",
                    annot.ext=annot,
                    allowMultiOverlap = TRUE,
                    isPairedEnd=TRUE,
                    requireBothEndsMapped=TRUE,
                    countMultiMappingReads=FALSE,
                    minOverlap = 10,
                    nthreads=8)
fc1$counts[rownames(fc1$counts) == 'orf14869',]

gtf <- file.path('~/RNASeq/YBR/scer.gtf')
fc <- featureCounts(files="~/RNASeq/YBR/fastq_files/YBR3/rawdata/quant_star/3SP1_quantAligned.sortedByCoord.out.bam",
                    annot.ext=gtf,isGTFAnnotationFile=TRUE,
                    isPairedEnd=TRUE,
                    requireBothEndsMapped=TRUE,
                    countMultiMappingReads=FALSE,
                    minOverlap = 10,
                    nthreads=8)
fc$counts[rownames(fc$counts) == 'YBR196C-A',]


hi <- merge(fc1$annotation, fc$annotation, by = c('Chr','Start','End','Strand','Length'))
hello <- fc$annotation[fc$annotation$GeneID %notin% hi$GeneID.y,]
hello2 <- hello[!str_detect(hello$Chr, ';'),-6]
hello3 <- hello[str_detect(hello$Chr, ';'),-6]

hello3$Chr <- str_split(hello3$Chr, ';', simplify = T)[,1]
hello3$Start <- str_split(hello3$Start, ';', simplify = T)[,1]
hello3$End <- str_split(hello3$End, ';', simplify = T)[,1]
hello3$Strand <- str_split(hello3$Strand, ';', simplify = T)[,1]


annot2 <- rbind(annot,hello2,hello3)
annot2 <- merge(annot2, fc$annotation[,-6], by = c('Chr','Start','End','Strand'), suffix = c('','_SGD'), all.x = T)
annot2 <- annot2[,c(5,1,2,3,4,6)]

write.csv(annot2, file = '~/RNASeq/YBR/translatome_annotation.csv')

