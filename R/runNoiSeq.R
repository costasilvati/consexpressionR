#' Execute NOISeq gene Expression analysis
#'
#' @param countMatrix either a matrix of raw (read) counts.
#' @param designExperiment replicate and treatment by samples
#' @param kParam Counts equal to 0 are replaced by k. By default, k = 0.5
#' @param factorParam A string indicating the name of factor whose levels are the conditions to be compared.
#' @param lcParam Length correction is done by dividing expression by length^lc. By default, lc = 0
#' @param replicatesParam In this argument, the type of replicates to be used is defined: "technical", "biological" or "no" replicates. By default, "technical" replicates option is chosen.
#' @param normParm Normalization method t can be one of "rpkm" (default), "uqua" (upper quartile), "tmm" (trimmed mean of M) or "n" (no normalization).
#'
#' @return A data.frame with output of noiseqOut@results
#' @export
#'
#' @examples
#' groupNameModel = c("BM","JJ")
#' numberReplicsModel = 3
#' designExperimentModel <- rep(groupNameModel, each = numberReplicsModel)
#' toolResult <- NULL
#' toolResult$noiseq <- runNoiSeq(countMatrix = gse95077,
#'                                designExperiment = designExperimentModel)
runNoiSeq <- function (countMatrix,
                       designExperiment,
                       normParm = "rpkm",
                       kParam = 0.5,
                       factorParam="Tissue",
                       lcParam = 0,
                       replicatesParam = "technical"){
  myfactors = data.frame(Tissue=c(designExperiment))
  print(myfactors$Tissue)
  mydata <- NOISeq::readData(data=countMatrix, factors=myfactors)
  print(utils::head(mydata$Tissue))
  #pnr?? nss?? v??
  mynoiseq = NOISeq::noiseq(mydata,
                            norm = normParm,
                            k= kParam,
                            factor=factorParam,
                            lc= lcParam,
                            replicates=replicatesParam)
  result <- mynoiseq@results[[1]]
  return(result)
}
