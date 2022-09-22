library(gtools)
library(dplyr)
library(stringi)
library(stringr)
library(openxlsx)

lof.mutants <- read.csv("~/R/Projects/ybr/input/ybr_lof_mutants.csv", header = T, stringsAsFactors = F)

##### THOUGHT EXPT 1
thought_experiment <- function(mutants) {
  lof.mutants <- read.csv("~/R/Projects/ybr/input/ybr_lof_mutants.csv", header = T, stringsAsFactors = F)
  
  temp <- lof.mutants[lof.mutants$Mutant %in% mutants,]
  cause <- data.frame(as.list(colMeans(temp[,2:11])))
  pheno.same <- colnames(cause)[cause[1,] == 1]
  pheno.diff <- colnames(cause)[cause[1,] > 0 & cause[1,] < 1]
  
  print(sprintf('%s ::    Same Phenotype   :: %s', paste(mutants, collapse = " | "), paste(pheno.same, collapse = " | ")))
  print(sprintf('%s :: Different Phenotype :: %s', paste(mutants, collapse = " | "), paste(pheno.diff, collapse = " | ")))
} 

thought_experiment(c('ATG1','ATG12'))

##### THOUGHT EXPT 2
phenotypes <- c('neutral', 'deleterious','beneficial')

mutant_options <- unique(lof.mutants$Mutant[lof.mutants$Mutant != 'WT'])

mutant_list <- list(c('ATG12'),
  c('ATG12','STOP23','STOP23C'),
  c('ATG12','STOP23','STOP23C','DEL'),
  c('ATG12','DEL'),
  c('ATG1', 'ATG12', 'DEL', 'STOP23', 'STOP23C'),
  c('ATG1', 'ATG12', 'STOP23', 'STOP23C'),
  c('ATG12','STOP23','STOP23C','STOP2','STOP2C')
  )

all_results <- list()
success_rate <- NULL
ii <- 1
# for (n_mutants in seq(2,13)) {
#   mutant_list <- permutations(n=13,r=n_mutants,
#                               v=mutant_options,
#                               repeats.allowed=F)
  # for (l in seq(1,dim(mutant_list)[1])) {
    # mutants <- mutant_list[l,]
for (mutants in mutant_list) {
    results <- permutations(n=length(phenotypes),r=length(mutants),
                            v=phenotypes,
                            repeats.allowed=T) %>%
      data.frame(stringsAsFactors = T)
    colnames(results) <- mutants
    # head(results)
    
    for (i in 1:dim(results)[1]) {
      temp_res <- results[i,]
      temp_del <- colnames(temp_res)[temp_res == 'deleterious']
      temp_neu <- colnames(temp_res)[temp_res == 'neutral']
      temp_ben <- colnames(temp_res)[temp_res == 'beneficial']
      
      temp <- lof.mutants[lof.mutants$Mutant %in% temp_del,]
      if (dim(temp)[1] > 0) {
        cause_del <- data.frame(as.list(colMeans(temp[,2:11])))
        temp <- lof.mutants[lof.mutants$Mutant %in% temp_neu,]
        cause_neu <- data.frame(as.list(colMeans(temp[,2:11])))
        temp <- lof.mutants[lof.mutants$Mutant %in% temp_ben,]
        
        if (dim(temp)[1] > 0) {
          cause_ben <- data.frame(as.list(colMeans(temp[,2:11])))
          pheno_del <- colnames(cause_del)[cause_del[1,] == 1]
          pheno_neu <- colnames(cause_neu)[cause_neu[1,] > 0]
          pheno_ben <- colnames(cause_ben)[cause_ben[1,] > 0]
          cause <- pheno_del[!(pheno_del %in% c(pheno_neu, pheno_ben))]
        } else {
          pheno_del <- colnames(cause_del)[cause_del[1,] == 1]
          pheno_neu <- colnames(cause_neu)[cause_neu[1,] > 0]
          cause <- pheno_del[!(pheno_del %in% pheno_neu)]
        }
        
        if (length(cause) > 0) {
          results[i,'Cause_of_LOF'] <- paste(cause, collapse = " | ")
        } else {
          results[i,'Cause_of_LOF'] <- 'inconclusive'
        }
      } else {
        results[i,'Cause_of_LOF'] <- 'inconclusive'
      }
      
      temp <- lof.mutants[lof.mutants$Mutant %in% temp_ben,]
      if (dim(temp)[1] > 0) {
        cause_ben <- data.frame(as.list(colMeans(temp[,2:11])))
        temp <- lof.mutants[lof.mutants$Mutant %in% temp_neu,]
        cause_neu <- data.frame(as.list(colMeans(temp[,2:11])))
        
        temp <- lof.mutants[lof.mutants$Mutant %in% temp_del,]
        if (dim(temp)[1] > 0) {
          cause_del <- data.frame(as.list(colMeans(temp[,2:11])))
          pheno_del <- colnames(cause_del)[cause_del[1,] > 0]
          pheno_neu <- colnames(cause_neu)[cause_neu[1,] > 0]
          pheno_ben <- colnames(cause_ben)[cause_ben[1,] == 1]
          cause <- pheno_ben[!(pheno_ben %in% c(pheno_neu, pheno_del))]
        } else {
          pheno_ben <- colnames(cause_del)[cause_ben[1,] == 1]
          pheno_neu <- colnames(cause_neu)[cause_neu[1,] > 0]
          cause <- pheno_ben[!(pheno_ben %in% pheno_neu)]
        }
        
        if (length(cause) > 0) {
          results[i,'Cause_of_GOF'] <- paste(cause, collapse = " | ")
        } else {
          results[i,'Cause_of_GOF'] <- 'inconclusive'
        }
      } else {
        results[i,'Cause_of_GOF'] <- 'inconclusive'
      }
    }
    success_rate$mutants[ii] <- paste(mutants, collapse = " | ")
    success_rate$lof_certainlyybr[ii] <- dim(results[results$Cause_of_LOF %in% c('YBR_Protein','YBR_Translation'),])[1]/
      dim(results)[1] * 100
    success_rate$lof_maybeybr[ii] <- dim(results[str_detect(results$Cause_of_LOF, 'YBR'),])[1]/
      dim(results)[1] * 100
    success_rate$lof_inconclusive[ii] <- dim(results[results$Cause_of_LOF == 'inconclusive',])[1]/
      dim(results)[1] * 100
    success_rate$gof_certainlyybr[ii] <- dim(results[results$Cause_of_GOF %in% c('YBR_Protein','YBR_Translation'),])[1]/
      dim(results)[1] * 100
    success_rate$gof_maybeybr[ii] <- dim(results[str_detect(results$Cause_of_GOF, 'YBR'),])[1]/
      dim(results)[1] * 100
    success_rate$gof_inconclusive[ii] <- dim(results[results$Cause_of_GOF == 'inconclusive',])[1]/
      dim(results)[1] * 100
    ii <- ii + 1
    
    all_results[paste(mutants, collapse = "_")] <- list(results)
  }
# }
success_rate <- data.frame(success_rate)

write.xlsx(success_rate, file = "output/success_rate.xlsx")
write.xlsx(all_results, file = "output/thought_experiment.xlsx")

