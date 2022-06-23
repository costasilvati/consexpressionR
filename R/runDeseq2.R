runDeseq2 <- function(tableCountPath, countMatrix,
                     designExperiment,
                     outFile="consexpression2_deseq2.csv",
                     sepCharacter="\t",
                     fitTypeDispersion="local"){

  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = coldata,
                                design= ~ batch + condition)
  dds <- DESeq(dds)
  resultsNames(dds) # lists the coefficients
  res <- results(dds, name="condition_trt_vs_untrt")
  # or to shrink log fold changes association with condition:
  res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")
}
