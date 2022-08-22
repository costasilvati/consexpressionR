runDeseq2 <- function(countMatrix,
                     designExperiment){
  dds <-DESeq2::DESeqDataSetFromMatrix(countData = countMatrix,
                                       colData = colnames(countMatrix),
                                       design= designExperiment)
  dds <- DESeq2::DESeq(dds)
  resultsNames(dds) # lists the coefficients
  res <- DESeq2::results(dds,
                         name=nameResults)
  return(res)
}
