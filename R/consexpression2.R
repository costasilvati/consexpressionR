#' Main function. Make expression analysis of 7 tools and print output by consenso
#'
#' @param numberReplics number of replicate (technical or biological) by sample
#' @param groupName text, name of samples or treatment
#' @param tableCountPath path to file csv that contains cout data or abundance data (local)
#' @param sepCharacter character used to split csv data, can be comma or tab
#' @param experimentName text, name of experiment
#' @param outDirPath path to write output, need be a directory (local)
#' @param methodNorm normalization method to be used (limma)
#' @param methodAdjPvalue correction method, a character string. Can be abbreviated (limma)
#' @param numberTopTable maximum number of genes to list (limma)
#' @param printResults logical variable: TRUE print report by each tool, FALSE print only consensus result
#'
#' @return
#' @export
#'
#' @examples
#' consexpression2(numberReplics=3,
#                 groupName= c("0b", "1b"),
#                 tableCountPath="tablecount.csv",
#                 sepCharacter="\t",
#                 experimentName="genericExperiment",
#                 outDirPath="consexpression2_results/",
#                 methodNorm = "TMM",
#                 methodAdjPvalue = "BH",
#                 numberTopTable = 1000000,
#                 printResults=FALSE)
consexpression2 <- function (numberReplics,
                             groupName,
                             tableCountPath,
                             sepCharacter="\t",
                             experimentName="genericExperiment",
                             outDirPath="consexpression2_results/",
                             methodNorm = "TMM",
                             methodAdjPvalue = "BH",
                             numberTopTable = 1000000,
                             printResults=FALSE){
    # validações?
    countMatrix <- readCountFile(tableCountPath,
                                 sepCharacter)
    designExperiment <- rep(groupName,
                            each = numberReplics)
    result <- NULL
    #result$bayseq<-runBaySeq(countMatrix,
     #                        designExperiment)

    result$edger<-runEdger(countMatrix,
                          numberReplics,
                          designExperiment)

    result$limma<-runLimma(countMatrix,
                          numberReplics,
                          designExperiment,
                          methodNorm,
                          methodAdjPvalue,
                          numberTopTable)

    result$noiseq<-runNoiSeq(countMatrix,
                            designExperiment)

    result$ebseq <- runEbseq(countMatrix,
                             designExperiment)
    # DESeq2 ??
    result$desq2 <- runDeseq2(countMatrix,
                              designExperiment)
    # SAMSeq só para count data
    result$samseq<-runSamSeq(countMatrix,
                             numberReplics,
                             designExperiment)
}
