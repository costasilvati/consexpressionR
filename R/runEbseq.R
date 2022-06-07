#' Execute EBSeq gene Expression analysis
#'
#'@param countMatrix either a matrix of raw (read) counts.
#' @param designExperiment replicate and treatment by samples
#' @param outFile path to write EBSeq report in csv file format
#' @param fdr
#' @param outFile
#' @param sepCharacter
#'
#' @return EBSeq report in data Frame fromat
#' @export
#'
#' @examples
runEbseq <- function(countMatrix, designExperiment, fdr=0.05,outFile="consexpression2_noiseq.csv", sepCharacter){
  Sizes<-MedianNorm(countMatrix)
  EBOut=EBTest(Data=m,
               Conditions<-as.factor(rep(c(grup),each=str(self._replic))),
               sizeFactors=Sizes, maxround=5)
  EBDERes=GetDEResults(EBOut, FDR=fdr)
  write.table(EBDERes$Status, file=outFile, sep=sepCharacter, quote=FALSE)
}

