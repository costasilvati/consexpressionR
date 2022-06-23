#' Execute EBSeq gene Expression analysis
#'
#'@param countMatrix either a matrix of raw (read) counts.
#' @param designExperiment replicate and treatment by samples
#' @param outFile path to write EBSeq report in csv file format
#' @param fdr False Discovery Rate cutt off
#' @param outFile path to write EBSeq report in csv file format
#' @param sepCharacter pattern to seprare file csv
#' @param maxRound
#'
#' @return EBSeq report in data Frame fromat
#' @export
#'
#' @examples
runEbseq <- function(countMatrix, designExperiment, maxRound=5, fdr=0.05,outFile="consexpression2_ebseq.csv", sepCharacter){
  sizes<-MedianNorm(countMatrix)
  ebOut<-EBTest(Data=countMatrix,
               Conditions=as.factor(designExperiment),
               sizeFactors=sizes, maxround=maxRound)
  results<-GetDEResults(ebOut, FDR=fdr)
  write.table(results$Status, file=outFile, sep=sepCharacter, quote=FALSE)
  return(results$Status)
}

