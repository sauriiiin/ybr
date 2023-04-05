##### INITIALIZE
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 10/19/2022

#### PACKAGES
library(sva)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(edgeR)
library(UpSetR)
library(dplyr)
library(stringr)
library(genefilter)
library(apeglm)
library(clusterProfiler)
library(GO.db)
library(org.Sc.sgd.db)
library(reshape2)
library(RMariaDB)
library(DOSE)
library(enrichplot)
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")

conn <- initialize.sql("saurin_test")

`%notin%` <- Negate(`%in%`)

fig.path <- 'output/figs/'

##### FUNCTIONS
run_edgeR = function(data, group_a_name, group_a_samples, group_b_samples, group_b_name){
  #create a list of all samples for this current comparison
  samples_for_comparison = c(group_a_samples,group_b_samples)
  
  #define the class factor for this pair of sample sets
  class = factor(c(rep(group_a_name,length(group_a_samples)),rep(group_b_name,length(group_b_samples))))
  
  #create a simplified data matrix for only these samples
  rawdata = data[,samples_for_comparison]
  
  #store gene names for later
  genes = rownames(data)
  
  #make DGElist object
  y = DGEList(counts=rawdata, genes=genes, group=class)
  
  #perform TMM normalization
  y <- calcNormFactors(y)
  
  #estimate dispersion & differential expression test
  if (length(samples_for_comparison) > 2) {
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y)
    et <- exactTest(y, pair = c(group_a_name, group_b_name))
  } else {
    et <- exactTest(y, pair = c(group_a_name, group_b_name), dispersion = 0.2^2)
  }
  
  #print number of up/down significant genes at FDR = 0.05 significance level and store the DE status in a new variable (de)
  summary(de <- decideTestsDGE(et, p.value=.05))
  summary(de <- decideTestsDGE(et, adjust.method="BH", p.value=.05))
  
  #create a matrix of significant DE genes
  mat <- data.frame(
    orf_name = genes,
    log10p = log10(et$table$PValue),
    fdr = p.adjust(et$table$PValue, method = 'BH'),
    lfc = et$table$logFC,
    contrast = sprintf('%s_vs_%s', group_b_name, group_a_name)
  )
  
  #order by log fold change
  mat <- mat[order(mat$lfc, decreasing = T),]
  row.names(mat) <- NULL
  
  #fix the issue where corrected p-values that are 0 become -Inf upon log10 conversion
  x <- mat$log10p
  lowest_pvalue <- min(x[which(!x == -Inf)])
  i = which(mat[,"log10p"] == -Inf)
  mat[i,"log10p"] = lowest_pvalue
  
  return(mat)
}

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

