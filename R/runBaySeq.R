#' Execute baySeq gene Expression analysis
#'
#' @param countMatrix either a matrix of raw (read) counts.
#' @param designExperiment replicate and treatment by samples
#' @param out path to write baySeq report in csv file format
#'
#' @return baySeq report in data Frame fromat
#' @export
#' @import baySeq
#'
#' @examples
#' bayseqResult<-runBaySeq(countMatrix, replicates, createNameFileOutput(outDirPath,experimentName,execName='baySeq'))
runBaySeq <- function  (countMatrix,
                        groupName,
                        numberReplics,
                        sampleSize = 1000,
                        estimationMethod = "QL",
                        clusters=2,
                        equalDispersion = TRUE,
                        priorProbabilities=c(0.5, 0.5), #0.5 por grupo??
                        reestimationType = "BIC",
                        topCountGroup = "DE",
                        numberTopCount=65000){
    if(require("parallel")) cl <- parallel::makeCluster(clusters) else cl <- NULL
    # rep("experiment", length(designExperiment)
    groups <- list(NDE = colnames(countMatrix),
                   DE = rep(groupName, each=numberReplics))
    CD <- new("countData",
              data = as.matrix(countMatrix),
              replicates = rep(groupName, each=numberReplics),
              groups = groups)
    baySeq::libsizes(CD) <- baySeq::getLibsizes(CD)
    CD <- baySeq::getPriors.NB(CD,
                               samplesize = sampleSize,
                               estimation = estimationMethod,
                               cl = cl,
                               equalDispersions = equalDispersion)
    CD <- baySeq::getLikelihoods(CD,
                                 prs=priorProbabilities,
                                 pET=reestimationType,
                                 cl=cl)
    result <- baySeq::topCounts(CD,
                        group = topCountGroup,
                        number = numberTopCount)
    return (result)
}
