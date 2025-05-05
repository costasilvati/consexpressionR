#' Execute NOISeq gene Expression analysis
#'
#' @param countMatrix either a matrix of raw (read) counts.
#' @param designExperiment replicate and treatment by samples.
#' @param kParam Counts equal to 0 are replaced by k. By default, k = 0.5
#' @param factorParam A string indicating the name of factor whose levels are the conditions to be compared.
#' @param lcParam Length correction is done by dividing expression by length^lc. By default, lc = 0
#' @param replicatesParam In this argument, the type of replicates to be used is defined: "technical", "biological" or "no" replicates. By default, "technical" replicates option is chosen.
#' @param normParm Normalization method t can be one of "rpkm" (default), "uqua" (upper quartile), "tmm" (trimmed mean of M) or "n" (no normalization).
#' @param condExp A vector containing the two conditions to be compared by the differential expression algorithm (needed when the factor contains more than 2 different conditions).
#' @param groups text list separated by comma, name of samples or treatment.
#'
#' @return A data.frame with output of noiseqOut@results
#' @export
#'
#' @examples
#' data(gse95077)
#' treats = c("BM", "JJ")
#' toolResult <- NULL
#' toolResult$noiseq <- runNoiSeq(countMatrix = gse95077,
#'                                designExperiment = rep(treats, each = 3))
runNoiSeq <- function (countMatrix,
                       designExperiment,
                       groups = c(""),
                       normParm = "rpkm",
                       kParam = 0.5,
                       factorParam=c("Tissue"),
                       lcParam = 0,
                       replicatesParam = "technical",
                       condExp = c("")){

  myfactors = data.frame(Tissue=c(designExperiment))
  mydata <- NOISeq::readData(data=countMatrix, factors=myfactors)
  if(length(groups) > 2 && length(condExp) == 2){
    if(all(condExp %in% groups)){
      mynoiseq = NOISeq::noiseq(mydata,
                                norm = normParm,
                                k= kParam,
                                factor=factorParam,
                                lc= lcParam,
                                replicates=replicatesParam,
                                conditions = condExp)
    }else{
      print("Values in conditions need be founded in groups")
      return( NULL)
    }
  }else{
    mynoiseq = NOISeq::noiseq(mydata,
                              norm = normParm,
                              k= kParam,
                              factor=factorParam,
                              lc= lcParam,
                              replicates=replicatesParam)
  }
  #pnr?? nss?? v??
  result <- mynoiseq@results[[1]]
  return(result)
}
