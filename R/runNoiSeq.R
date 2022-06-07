#' Execute NOISeq gene Expression analysis
#'
#'
#' @param countMatrix either a matrix of raw (read) counts.
#' @param designExperiment replicate and treatment by samples
#' @param noiseqOutPtah path to write EBSeq report in csv file format
#'
#' @return NOISeq report in data Frame fromat
#' @export
#'
#' @examples
#' noiseqResult<-runNOISeq(countMatrix, designExperiment, createNameFileOutput(outDirPath,experimentName,execName='NOISeq'))
runNoiSeq <- function (countMatrix, designExperiment, noiseqOutPtah){
    myfactors = data.frame(Tissue=c(designExperiment))
    mydata <- readData(data=countMatrix, factors=myfactors)
    # k?? lc?? factor??
    mynoiseq = noiseq(mydata, k=0.5, factor="Tissue", lc=1, replicates="technical")
    result <- mynoiseq@results[[1]]
    write.csv(result, file=noiseqOutPtah, sep="\t", quote=FALSE)
    return(mynoiseq)
}
