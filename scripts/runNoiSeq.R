runNoiSeq <- function (countMatrix, designExperiment, noiseqOutPtah){
    library("NOISeq")
    myfactors = data.frame(Tissue=c(designExperiment))
    mydata <- readData(data=countMatrix, factors=myfactors)
    # k?? lc?? factor??
    mynoiseq = noiseq(mydata, k=0.5, factor="Tissue", lc=1, replicates="technical")
    result <- mynoiseq@results[[1]]
    write.csv(result, file=noiseqOutPtah, sep="\t", quote=FALSE)
    return(mynoiseq)
}