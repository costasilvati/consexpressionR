#' Execute EBSeq gene Expression analysis
#'
#' @param countMatrix either a matrix of raw (read) counts.
#' @param designExperiment replicate and treatment by samples
#' @param fdr parameter used in EBTest function: fdr False Discovery Rate cutt off
#' @param maxRound parameter used in EBTest function: Number of iterations. The default value is 50.
#' @param methodDeResults parameter used in GetDEResults function: "robust" or "classic". Using the "robust" option, EBSeq is more robust to genes with outliers and genes with extremely small variances. Using the "classic" option, the results will be more comparable to those obtained by using the GetPPMat() function from earlier version (<= 1.7.0) of EBSeq
#' @param groups text, name of samples or treatment
#' @param ppThreshold posterior Probability Threshold
#' @return EBSeq report in data Frame fromat
#' @export
#'
#' @examples
#' data(gse95077)
#' treats = c("BM", "JJ")
#' designExperimentModel <- rep(treats, each = 3)
#' toolResult <- NULL
#' toolResult$ebseq <- runEbseq(countMatrix = gse95077,designExperiment = designExperimentModel,groups = treats)
runEbseq <- function(countMatrix,designExperiment,fdr=0.05,ppThreshold = 0.8,maxRound = 50,
                     methodDeResults = "robust",groups = c("")){
  sizes<- EBSeq::MedianNorm(countMatrix)
  if(length(groups) > 2){
    cond <- as.factor(designExperiment)
    ebseqOutMultiTest <- EBSeq::EBMultiTest(Data = as.matrix(countMatrix),Conditions = cond,sizeFactors = sizes,maxround = maxRound)
    ebseqOutMultiPP <- EBSeq::GetMultiPP(ebseqOutMultiTest)
    linhas_selecionadas <- apply(ebseqOutMultiPP$PP, 1, function(row) any(row >= ppThreshold))
    ebseqOutMultiPP$difExp <- ifelse(linhas_selecionadas, "DE", "EE")
    result <- as.data.frame(ebseqOutMultiPP$difExp)
    colnames(result)[ncol(result)] <- "DEFound"
    return(result)
  }else{
    ebOut<-EBSeq::EBTest(Data=as.matrix(countMatrix),Conditions=as.factor(designExperiment),sizeFactors=sizes, maxround = maxRound)
    ebOutResult <- EBSeq::GetDEResults(ebOut,FDR=fdr)
    result <- as.data.frame(ebOutResult$Status)
    colnames(result)[ncol(result)] <- "DEFound"
    return(result)
  }
}

