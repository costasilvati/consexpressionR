runBaySeq <- function  (countMatrix, designExperiment, out){
    library("baySeq")
    if(require("parallel")) cl <- makeCluster(4) else cl <- NULL
    groups <- list(NDE = designExperiment, DE = designExperiment)
    CD <- new("countData", data = countMatrix, replicates = designExperiment, groups = groups)
    libsizes(CD) <- getLibsizes(CD)
    CD <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl, equalDispersions = TRUE)
    CD <- getLikelihoods(CD, prs=c(0.5, 0.5), pET="BIC", cl=cl)
    write.table(topCounts(CD, group = "DE", number = 65000, normaliseData = TRUE), out, sep="\t", quote = FALSE)
    return (CD)
}