# Title     : TabularCounts
# Objective : Tabular count or quant results of mapping reads
# Created by: julianacostasilva
# Created on: 01/11/21
local_path <- '/Users/julianacostasilva/Library/CloudStorage/OneDrive-Pessoal/ProjetoDoutorado/bioconvergencia/reads_RNApa/kallisto_quant_align_apa_1B_0B/gonoda/'
pattern_find <- "*.tsv"
file_out <- paste(local_path,"gonoda_tpm.csv", sep = "")
column_tpm <- 'target_id'
setwd(local_path)

list_dir <- list.files(path = ".", pattern =pattern_find, all.files = FALSE,
           full.names = FALSE, recursive = TRUE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
flag <- TRUE

for(sample in list_dir){
  data <- data.frame(read.csv(sample, header = 1, sep = '\t', row.names = column_tpm))
  coluname <- substring(sample, 1, 5)
  #renomeia tpm para amostra
  colnames(data)[4] <- coluname
  # remove colunas
  data <- data[-c(1:3)]

  if (flag){
    new_tabular <- data
    flag <- FALSE
  }else{
    new_tabular <- cbind.data.frame(new_tabular, data)
    write.csv(new_tabular, file_out, strin)
    write
  }
}
