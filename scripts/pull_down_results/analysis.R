
ms_overlap <- function(combine_replicates, test_group, control_group) {

  # combine_replicates <- 'N'
  # test_group <- c("attempt_2_WT_YBR_HA") #"attempt_2_WT_YBR_mNG", "attempt_2_doa10D_YBR_mNG", "attempt_1_doa10D_YBR_mNG",
  # control_group <- c("attempt_2_WT_HA") #"attempt_2_WT_mNG", "attempt_2_doa10D_mNG", "attempt_1_WT_mNG",
  # 
  library(dplyr)
  library(stringr)
  library(reshape2)
  library(UpSetR)
  library(clusterProfiler)
  library(GO.db)
  library(org.Sc.sgd.db)
  source('/home/sbp29/R/Projects/ybr/scripts/run_edgeR.R')
  
  `%notin%` <- Negate(`%in%`)
  all_genes <- read.csv(file = '/home/sbp29/R/Projects/ybr/input/ms_data/scer_gene_IDs.csv')
  
  
  pdone <- read.csv(file = '/home/sbp29/R/Projects/ybr/scripts/pull_down_results/pull_down_one.csv')
  pdone.cnts <- read.csv(file = '/home/sbp29/R/Projects/ybr/scripts/pull_down_results/pull_down_one_cnts.csv')
  pdtwo <- read.csv(file = '/home/sbp29/R/Projects/ybr/scripts/pull_down_results/pull_down_two.csv')
  pdtwo.cnts <- read.csv(file = '/home/sbp29/R/Projects/ybr/scripts/pull_down_results/pull_down_two_cnts.csv')
  sample_info <- read.csv(file = '/home/sbp29/R/Projects/ybr/input/ms_data/221225_sample_info.csv')
  sample_info$Sample.ID.number <- str_replace_all(sample_info$Sample.ID.number, '-', '.')
  
  pddat <- rbind(melt(pdone, id.vars = 'Gene', variable.name = 'sample', value.name = 'protein_present'),
                 melt(pdtwo, id.vars = 'Gene', variable.name = 'sample', value.name = 'protein_present'))
  pddat <- merge(pddat, sample_info[,c(1,2,8:10)], by.x = 'sample', by.y = 'Sample.ID.number')
  pddat <- pddat %>%
    mutate(ID = paste('attempt', Attempt, background, ybr_tag, sep = '_'))
  
  pddat.cnts <- rbind(melt(pdone.cnts, id.vars = 'Gene', variable.name = 'sample', value.name = 'peptide_cnts'),
                 melt(pdtwo.cnts, id.vars = 'Gene', variable.name = 'sample', value.name = 'peptide_cnts'))
  pddat.cnts <- merge(pddat.cnts, sample_info[,c(1,2,8:10)], by.x = 'sample', by.y = 'Sample.ID.number')
  pddat.cnts <- pddat.cnts %>%
    mutate(ID = paste('attempt', Attempt, background, ybr_tag, sep = '_'))
  pddat.cnts <- pddat.cnts %>%
    group_by(Attempt, background, ybr_tag, replicate, ID, Gene) %>%
    summarise(peptide_cnts = sum(peptide_cnts), .groups = 'keep') %>%
    filter(peptide_cnts > 0) %>%
    data.frame()
  
  
  ##### UPSET PLOT
  temp.upset <- data.frame(Gene = unique(pddat$Gene))
  if (combine_replicates == 'Y') {
    for (g in c(test_group, control_group)) {
      temp.group <- pddat %>%
        filter(ID %in% g) %>%
        group_by(Gene) %>%
        summarise(protein_present = sum(protein_present)) %>%
        mutate(protein_present = case_when(protein_present > 0 ~ 1,
                                           protein_present == 0 ~ 0)) %>%
        data.frame()
      
      temp.upset <- merge(temp.upset, temp.group, by = 'Gene', all.x = T)
    }
    colnames(temp.upset) <- c('Gene', c(test_group, control_group))
    temp.upset[is.na(temp.upset)] <- 0
    upset_plot <- upset(temp.upset, nsets = ncol(temp.upset) - 1, order.by = 'freq', text.scale = 2)
  } else {
    col.names <- NULL
    for (g in c(test_group, control_group)) {
      for (r in unique(pddat$replicate[pddat$ID == g])) {
        temp.group <- pddat %>%
          filter(ID %in% g, replicate %in% r) %>%
          group_by(Gene) %>%
          summarise(protein_present = sum(protein_present)) %>%
          mutate(protein_present = case_when(protein_present > 0 ~ 1,
                                             protein_present == 0 ~ 0)) %>%
          data.frame()
        temp.upset <- merge(temp.upset, temp.group, by = 'Gene', all.x = T)
        col.names <- c(col.names, paste(g,r,sep='_'))
      }
    }
    colnames(temp.upset) <- c('Gene', col.names)
    temp.upset[is.na(temp.upset)] <- 0
    upset_plot <- upset(temp.upset, nsets = ncol(temp.upset) - 1, order.by = 'freq', text.scale = 2)
  }
  
  ###### TEST GROUP ONLY GENES
  # ID count
  n_groups <- length(test_group)
  n_groups_reps <- pddat.cnts %>%
    filter(ID %in% test_group) %>%
    group_by(ID, replicate) %>% count() %>% nrow()
  
  # replicate count per ID
  test_group.reps <- pddat %>%
    filter(ID %in% test_group, protein_present == 1) %>%
    group_by(ID, replicate) %>%
    count() %>% 
    group_by(ID) %>%
    count() %>% data.frame()
  
  # genes in control group
  control_genes <- pddat %>%
    filter(ID %in% control_group, protein_present == 1) %>%
    group_by(ID, replicate, Gene) %>%
    count() %>% 
    group_by(ID, replicate, Gene) %>%
    count() %>% data.frame()
  control_genes <- unique(control_genes$Gene)
  
  control_genes.pep_cnts <- pddat.cnts %>%
    filter(ID %in% control_group) %>%
    group_by(Gene) %>%
    summarise(peptide_cnts_control = sum(peptide_cnts, na.rm = T)) %>%
    data.frame()
  
  # filter genes in control group from test group
  # genes per ID per replicate
  test_genes.perID.perRep <- pddat %>%
    filter(ID %in% test_group, protein_present == 1) %>%
    filter(Gene %notin% unique(control_genes)) %>%
    group_by(ID, replicate, Gene) %>%
    count() %>%
    group_by(ID, replicate, Gene) %>%
    count() %>% data.frame()
  
  test_genes.pep_cnts <- pddat.cnts %>%
    filter(ID %in% test_group) %>%
    group_by(Gene) %>%
    summarise(peptide_cnts_test = sum(peptide_cnts, na.rm = T)) %>%
    data.frame()

  # test group only genes per ID (replicate combined)
  test_genes.perID <- test_genes.perID.perRep %>%
    group_by(ID, Gene) %>%
    count() %>% data.frame()
  
  # test group only genes count per ID
  test_genes.perID.count <- test_genes.perID %>%
    group_by(ID) %>%
    count() %>% data.frame()
  
  # test group genes count per ID
  test_genes.perID.allcount <- pddat %>%
    filter(ID %in% test_group, protein_present == 1) %>%
    group_by(ID, replicate, Gene) %>%
    count() %>%
    group_by(ID, replicate, Gene) %>%
    count() %>% data.frame() %>%
    group_by(ID, Gene) %>%
    count() %>% 
    group_by(ID) %>%
    count() %>% data.frame()
  
  # test group only genes that overlap
  test_genes.overlap <- test_genes.perID %>% 
    group_by(ID, Gene) %>%
    count() %>% 
    group_by(Gene) %>%
    count() %>% 
    filter(n == n_groups) %>% data.frame()
  
  test_genes.overlap.total <- pddat.cnts %>%
    filter(ID %in% test_group) %>%
    group_by(Gene) %>% count() %>%
    filter(n == n_groups_reps) %>%
    data.frame()
  
  # avg. pep. count of test group only genes that overlap
  test_genes.overlap.cnts <- pddat.cnts %>%
    filter(ID %in% test_group, Gene %in% test_genes.overlap$Gene) %>%
    group_by(Gene) %>%
    summarise(total_peptide_cnts = sum(peptide_cnts, na.rm = T),
              avg_peptide_cnts = mean(peptide_cnts, na.rm = T), .groups = 'keep') %>%
    data.frame()
  
  colnames(test_group.reps) <- c('ID','reps')
  colnames(test_genes.perID.allcount) <- c('ID','total_proteins_found')
  colnames(test_genes.perID.count) <- c('ID','total_proteins_found_minus_control')
  
  overlap_results <- merge(test_group.reps, test_genes.perID.allcount, by = 'ID', all = T)
  overlap_results <- merge(overlap_results, test_genes.perID.count, by = 'ID', all = T)
  overlap_results$overlap <- nrow(test_genes.overlap)
  
  test_genes.perID$overlap[test_genes.perID$Gene %in% test_genes.overlap$Gene] <- 1
  test_genes.perID$overlap[is.na(test_genes.perID$overlap)] <- 0
  
  
  if (sum(test_group.reps$reps > 1) >= 1) {
    # test group only genes that overlap per replicate when replicate > 1
    test_genes.perID.overlap <- merge(test_genes.perID, test_group.reps,
                                      by = 'ID', suffixes = c('','_max')) %>%
      filter(reps > 1, n == reps)
    
    # test group only genes that overlap per replicate and then between IDs
    test_genes.overlap.overlap <- merge(test_genes.perID, test_group.reps,
                                        by = 'ID', suffixes = c('','_max')) %>%
      filter(n == reps) %>%
      group_by(ID, Gene) %>%
      count() %>% 
      group_by(Gene) %>%
      count() %>% 
      filter(n == n_groups) %>% data.frame()
    
    # sprintf('%s has/have more than 1 replicate.', paste(test_group.reps$ID[test_group.reps$n > 1], collapse = ', '))
    
    test_genes.perID.overlap.count <- test_genes.perID.overlap %>%
      group_by(ID) %>%
      count() %>% data.frame()
    colnames(test_genes.perID.overlap.count) <- c('ID', 'n_genes_overlap_within_replicates')
    
    overlap_results <- merge(overlap_results, test_genes.perID.overlap.count, by = 'ID', all = T)
    overlap_results$overlap_with_replicate_wise_overlap <- nrow(test_genes.overlap.overlap)
    
    test_genes.perID$overlap_with_replicate_wise_overlap[test_genes.perID$Gene %in% test_genes.overlap.overlap$Gene] <- 1
    test_genes.perID$overlap_with_replicate_wise_overlap[is.na(test_genes.perID$overlap_with_replicate_wise_overlap)] <- 0
  }
  
  ##### ENRICHMENT ANALYSIS
  # using edgeR
  
  control_test.pep_cnts <- merge(control_genes.pep_cnts, test_genes.pep_cnts, by = 'Gene', all = T)
  control_test.pep_cnts[is.na(control_test.pep_cnts)] <- 0
  
  edgeR_in <- control_test.pep_cnts
  row.names(edgeR_in) <- edgeR_in$Gene
  edgeR_in <- edgeR_in[-1]
  
  edgeR_out <- run_edgeR(data=edgeR_in,
                         group_a_name='Control', group_a_samples='peptide_cnts_control',
                         group_b_name='YBR', group_b_samples='peptide_cnts_test')
  edgeR_out <- merge(control_test.pep_cnts, edgeR_out, by.x = 'Gene', by.y = 'orf_name', all = T)
  
  edgeR_out$raw_lfc <- log2(edgeR_out$peptide_cnts_test/edgeR_out$peptide_cnts_control)
  edgeR_out$norm_lfc <- log2((edgeR_out$peptide_cnts_test/sum(edgeR_out$peptide_cnts_test))/
                               (edgeR_out$peptide_cnts_control/sum(edgeR_out$peptide_cnts_control)))
  # edgeR_out$test_only_overlap[edgeR_out$Gene %in% test_genes.overlap$Gene] <- 1
  # edgeR_out$test_only_overlap[is.na(edgeR_out$test_only_overlap)] <- 0
  edgeR_out$test_overlap[edgeR_out$Gene %in% test_genes.overlap.total$Gene] <- 1
  edgeR_out$test_overlap[is.na(edgeR_out$test_overlap)] <- 0
  edgeR_out <- edgeR_out[order(edgeR_out$lfc, edgeR_out$peptide_cnts_test, decreasing = T),]
  
  ##### GO KEGG ENRICHMENT FROM THE TEST GROUP
  # test_genes.overlap
  gene_list <- all_genes[(all_genes$ORF %in% test_genes.overlap$Gene) |
                           (all_genes$GENENAME %in% test_genes.overlap$Gene),]
  
  goe <- enrichGO(gene          = gene_list$ENSEMBL,
                  universe      = all_genes$ENSEMBL,
                  OrgDb         = org.Sc.sgd.db,
                  keyType       = "ENSEMBL",
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05)
  goe <- data.frame(goe)
  
  kegg <- enrichKEGG(gene         = gene_list$ENSEMBL,
                     universe     = all_genes$ENSEMBL,
                     organism     = 'sce',
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff  = 0.05s)
  kegg <- data.frame(kegg)
  
  goe$GeneRatio <- as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,1])/
    as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,2])
  goe$BgRatio <- as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,1])/
    as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,2])
  goe$GO <- paste0(goe$ONTOLOGY, '_', goe$Description)
  goe <- goe[order(goe$GeneRatio, goe$qvalue, decreasing = T),]
  
  kegg$GeneRatio <- as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,1])/
    as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,2])
  kegg$BgRatio <- as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,1])/
    as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,2])
  kegg <- kegg[order(kegg$GeneRatio, kegg$qvalue, decreasing = T),]
  
  # head(goe)
  # head(kegg)
  
##### RETURN DATA
  output <- list(upset_plot, overlap_results, test_genes.perID, test_genes.overlap.cnts, edgeR_out, goe, kegg)
  return(output)
}








