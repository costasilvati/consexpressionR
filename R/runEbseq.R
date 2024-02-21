#' Execute EBSeq gene Expression analysis
#'
#' @param countMatrix either a matrix of raw (read) counts.
#' @param designExperiment replicate and treatment by samples
#' @param fdr False Discovery Rate cutt off
#' @param sepCharacter pattern to seprare file csv
#'
#' @return EBSeq report in data Frame fromat
#' @export
#'
#' @examples
runEbseq <- function(countMatrix,
                     designExperiment,
                     fdr=0.05,
                     sepCharacter){
  results <- NULL
  sizes<- EBSeq::MedianNorm(countMatrix)
  ebOut<-EBSeq::EBTest(Data=countMatrix,
                       Conditions=as.factor(designExperiment),
                       sizeFactors=sizes, maxround = 5)
  results<-EBSeq::GetDEResults(ebOut,
                               FDR=fdr)
  return(results$Status)
}

