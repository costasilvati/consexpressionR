
runVariancePatition <- function(countMatrix){
  # filter genes by number of counts
  isexpr <- rowSums(cpm(countMatrix) > 0.1) >= 5

  # Standard usage of limma/voom
  dge <- DGEList(countMatrix[isexpr, ])
  dge <- calcNormFactors(dge)

  # make this vignette faster by analyzing a subset of genes
  dge <- dge[1:1000, ]

  # Specify parallel processing parameters
  # this is used implicitly by dream() to run in parallel
  param <- SnowParam(1, "SOCK", progressbar = TRUE)

  # The variable to be tested must be a fixed effect
  form <- ~ Disease + (1 | Individual)

  # estimate weights using linear mixed model of dream
  vobjDream <- voomWithDreamWeights(dge, form, metadata, BPPARAM = param)

  # Fit the dream model on each gene
  # For the hypothesis testing, by default,
  # dream() uses the KR method for <= 20 samples,
  # otherwise it uses the Satterthwaite approximation
  fitmm <- dream(vobjDream, form, metadata)
  fitmm <- eBayes(fitmm)

}
