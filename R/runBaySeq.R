#' Execute baySeq gene Expression analysis
#'
#' @param countMatrix either a matrix of raw (read) counts.
#' @param designExperiment replicate and treatment by samples
#' @param out path to write baySeq report in csv file format
#'
#' @return baySeq report in data Frame fromat
#' @export
#'
#' @examples
#' bayseqResult<-runBaySeq(countMatrix, replicates, createNameFileOutput(outDirPath,experimentName,execName='baySeq'))
runBaySeq <- function  (countMatrix, designExperiment, out){
    if(require("parallel")) cl <- makeCluster(4) else cl <- NULL
    groups <- list(NDE = designExperiment, DE = designExperiment)
    CD <- new("countData", data = countMatrix, replicates = designExperiment, groups = groups)
    baySeq::libsizes(CD) <- baySeq::getLibsizes(CD)
    CD <- baySeq::getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl, equalDispersions = TRUE)
    CD <- baySeq::getLikelihoods(CD, prs=c(0.5, 0.5), pET="BIC", cl=cl)
    write.table(topCounts(CD, group = "DE", number = 65000, normaliseData = TRUE), out, sep="\t", quote = FALSE)
    return (CD)
}
