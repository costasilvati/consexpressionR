#' Execute EBSeq gene Expression analysis
#'
#' @param countMatrix either a matrix of raw (read) counts.
#' @param designExperiment replicate and treatment by samples
#' @param fdr parameter used in EBTest function: fdr False Discovery Rate cutt off
#' @param maxRound parameter used in EBTest function: Number of iterations. The default value is 50.
#' @param methodDeResults parameter used in GetDEResults function: "robust" or "classic". Using the "robust" option, EBSeq is more robust to genes with outliers and genes with extremely small variances. Using the "classic" option, the results will be more comparable to those obtained by using the GetPPMat() function from earlier version (<= 1.7.0) of EBSeq
#'
#' @return EBSeq report in data Frame fromat
#' @export
#'
#' @examples
#' groupNameModel = c("BM","JJ")
#' numberReplicsModel = 3
#' designExperimentModel <- rep(groupNameModel, each = numberReplicsModel)
#' toolResult <- NULL
#' toolResult$ebseq <- runEbseq(countMatrix = as.matrix(gse95077),
#'                               designExperiment = designExperimentModel)
runEbseq <- function(countMatrix,
                     designExperiment,
                     fdr=0.05,
                     maxRound = 50,
                     methodDeResults = "robust"){
  results <- NULL
  sizes<- EBSeq::MedianNorm(countMatrix)
  ebOut<-EBSeq::EBTest(Data=countMatrix,
                       Conditions=as.factor(designExperiment),
                       sizeFactors=sizes, maxround = maxRound)
  results<-EBSeq::GetDEResults(ebOut,
                               FDR=fdr)
  return(results$Status)
}

