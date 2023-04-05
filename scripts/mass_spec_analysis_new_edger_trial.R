
head(test_genes.pep_cnts)
head(control_genes.pep_cnts)

control_test.pep_cnts <- merge(control_genes.pep_cnts, test_genes.pep_cnts, by = 'Gene', all = T)
control_test.pep_cnts[is.na(control_test.pep_cnts)] <- 0

edgeR_in <- control_test.pep_cnts
row.names(edgeR_in) <- edgeR_in$Gene
edgeR_in <- edgeR_in[-1]
head(edgeR_in)

edgeR_out <- run_edgeR(data=edgeR_in,
          group_a_name='Control', group_a_samples='peptide_cnts_control',
          group_b_name='YBR', group_b_samples='peptide_cnts_test')
edgeR_out <- merge(control_test.pep_cnts, edgeR_out, by.x = 'Gene', by.y = 'orf_name', all = T)

edgeR_out$raw_lfc <- log2(edgeR_out$peptide_cnts_test/edgeR_out$peptide_cnts_control)
edgeR_out$norm_lfc <- log2((edgeR_out$peptide_cnts_test/sum(edgeR_out$peptide_cnts_test))/
                             (edgeR_out$peptide_cnts_control/sum(edgeR_out$peptide_cnts_control)))
edgeR_out$test_only_overlap[edgeR_out$Gene %in% test_genes.overlap$Gene] <- 1
edgeR_out$test_only_overlap[is.na(edgeR_out$test_only_overlap)] <- 0
edgeR_out$test_overlap[edgeR_out$Gene %in% test_genes.overlap.total$Gene] <- 1
edgeR_out$test_overlap[is.na(edgeR_out$test_overlap)] <- 0
head(edgeR_out)

# edgeR_out %>%
#   ggplot(aes(x = raw_lfc, y = norm_lfc)) +
#   geom_abline(col = 'red', linetype = 'dashed') +
#   geom_point() +
#   # geom_point(data = edgeR_out %>% filter(fdr <= 0.05),
#   #            aes(x = raw_lfc, y = norm_lfc),
#   #            col = 'red') +
#   coord_cartesian(xlim = c(-6,6),
#                   ylim = c(-6,6))
# 
# edgeR_out %>%
#   ggplot(aes(x = lfc, y = norm_lfc)) +
#   geom_abline(col = 'red', linetype = 'dashed') +
#   geom_point() +
#   geom_point(data = edgeR_out %>% filter(fdr <= 0.05),
#              aes(x = lfc, y = norm_lfc),
#              col = 'red') +
#   coord_cartesian(xlim = c(-10,10),
#                   ylim = c(-10,10))
# 
# edgeR_out %>%
#   ggplot(aes(x = peptide_cnts_test + peptide_cnts_control, y = lfc)) +
#   geom_point() + 
#   geom_point(data = edgeR_out %>% filter(fdr <= 0.05),
#              aes(x = peptide_cnts_test + peptide_cnts_control, y = lfc),
#              col = 'red') +
#   scale_x_log10()
# 
# edgeR_out %>%
#   ggplot(aes(x = peptide_cnts_test - peptide_cnts_control, y = lfc)) +
#   geom_point() + 
#   geom_point(data = edgeR_out %>% filter(fdr <= 0.05),
#              aes(x = peptide_cnts_test - peptide_cnts_control, y = lfc),
#              col = 'red')
# 
# edgeR_out %>%
#   ggplot(aes(x = lfc, y = -log10p)) +
#   geom_point() +
#   geom_point(data = edgeR_out %>% filter(fdr <= 0.05),
#              aes(x = lfc, y = -log10p),
#              col = 'red')
