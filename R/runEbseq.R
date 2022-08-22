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
runEbseq <- function(countMatrix,
                     designExperiment,
                     fdr=0.05,
                     sepCharacter){
  sizes<- EBSeq::MedianNorm(countMatrix)
  ebOut<-EBSeq::EBTest(Data=countMatrix,
                       Conditions=as.factor(designExperiment),
                       sizeFactors=sizes, maxround = 5)
  results<-EBSeq::GetDEResults(ebOut,
                               FDR=fdr)
  return(results$Status)
}

