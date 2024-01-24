#' Execute NOISeq gene Expression analysis
#'
#' @param countMatrix either a matrix of raw (read) counts.
#' @param designExperiment replicate and treatment by samples
#'
#' @return NOISeq report in data Frame fromat
#' @export
#'
#' @examples
#'
runNoiSeq <- function (countMatrix, designExperiment){
    myfactors = data.frame(Tissue=c(designExperiment))
    mydata <- NOISeq::readData(data=countMatrix, factors=myfactors)
    # k?? lc?? factor??
    mynoiseq = NOISeq::noiseq(mydata, k=0.5, factor="Tissue", lc=1, replicates="technical")
    result <- mynoiseq@results[[1]]
    return(result)
}
