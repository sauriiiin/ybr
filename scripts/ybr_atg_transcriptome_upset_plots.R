
#upset plots to look at overlap [Not that informative]
listInput1 = list("Batch1 WT vs Batch2 WT" = ybr1_wt_vs_ybr3_wt_uncorrected$orf_name[ybr1_wt_vs_ybr3_wt_uncorrected$fdr <= 0.05],
                  "Batch1 ATG vs Batch2 ATG" = ybr1_atg_vs_ybr3_atg_uncorrected$orf_name[ybr1_atg_vs_ybr3_atg_uncorrected$fdr <= 0.05])
listInput2 = list("Batch1 WT vs Batch2 WT" = ybr1_wt_vs_ybr3_wt_corrected$orf_name[ybr1_wt_vs_ybr3_wt_corrected$fdr <= 0.05],
                  "Batch1 ATG vs Batch2 ATG" = ybr1_atg_vs_ybr3_atg_corrected$orf_name[ybr1_atg_vs_ybr3_atg_corrected$fdr <= 0.05])
listInput3 = list("Batch1 WT vs Batch1 ATG" = ybr1_wt_vs_ybr1_atg_uncorrected$orf_name[ybr1_wt_vs_ybr1_atg_uncorrected$fdr <= 0.05],
                  "Batch2 WT vs Batch2 ATG" = ybr3_wt_vs_ybr3_atg_uncorrected$orf_name[ybr3_wt_vs_ybr3_atg_uncorrected$fdr <= 0.05])
listInput4 = list("Batch1 WT vs Batch1 ATG" = ybr1_wt_vs_ybr1_atg_corrected$orf_name[ybr1_wt_vs_ybr1_atg_corrected$fdr <= 0.05],
                  "Batch2 WT vs Batch2 ATG" = ybr3_wt_vs_ybr3_atg_corrected$orf_name[ybr3_wt_vs_ybr3_atg_corrected$fdr <= 0.05])
listInput5 = list("Batch1 WT vs Batch1 ATG" = ybr1_wt_vs_ybr1_atg_corrected$orf_name[ybr1_wt_vs_ybr1_atg_corrected$fdr <= 0.05],
                  "Batch2 WT vs Batch2 ATG" = ybr3_wt_vs_ybr3_atg_corrected$orf_name[ybr3_wt_vs_ybr3_atg_corrected$fdr <= 0.05],
                  "Combined WT vs Combined ATG" = ybr_wt_vs_ybr_atg_corrected$orf_name[ybr_wt_vs_ybr_atg_corrected$fdr <= 0.05])


# jpeg(filename = sprintf('%s/YBR_ATG_BATCH_COMPARE_UNCORRECTED.jpg',fig.path),
#      width = two.c, height = one.5c, res = 300, unit = 'mm', bg = 'white')
# upset(fromList(listInput1), order.by = "freq", point.size=3, set_size.scale_max = 6000, text.scale = 1)
# dev.off()
# 
# jpeg(filename = sprintf('%s/YBR_ATG_BATCH_COMPARE_CORRECTED.jpg',fig.path),
#      width = two.c, height = one.5c, res = 300, unit = 'mm', bg = 'white')
# upset(fromList(listInput2), order.by = "freq", point.size=3, set_size.scale_max = 6000, text.scale = 1)
# dev.off()
# 
# jpeg(filename = sprintf('%s/YBR_ATG_SAMPLE_COMPARE_UNCORRECTED.jpg',fig.path),
#      width = two.c, height = one.5c, res = 300, unit = 'mm', bg = 'white')
# upset(fromList(listInput3), order.by = "freq", point.size=3, set_size.scale_max = 6000, text.scale = 1)
# dev.off()
# 
# jpeg(filename = sprintf('%s/YBR_ATG_SAMPLE_COMPARE_CORRECTED.jpg',fig.path),
#      width = two.c, height = one.5c, res = 300, unit = 'mm', bg = 'white')
# upset(fromList(listInput4), order.by = "freq", point.size=3, set_size.scale_max = 6000, text.scale = 1)
# dev.off()
# 
# jpeg(filename = sprintf('%s/YBR_ATG_SAMPLE_COMPARE_ALL_CORRECTED.jpg',fig.path),
#      width = two.c, height = one.5c, res = 300, unit = 'mm', bg = 'white')
# upset(fromList(listInput5), order.by = "freq", point.size=5, set_size.scale_max = 6000, text.scale = 1)
# dev.off()