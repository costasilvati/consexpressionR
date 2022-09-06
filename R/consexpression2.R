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
#' @param kallistoReport
#' @param kallistoDir
#' @param kallistoSubDir
#' @param kallistoOut
#'
#' @return
#' @export
#'
#' @examples
#' consexpression2(numberReplics=3,
#                 groupName= c("X0b", "X1b"),
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
                             printResults=FALSE,
                             kallistoReport = "report.txt",
                             kallistoDir = "kallisto_quant",
                             kallistoSubDir = "expermient_kallisto",
                             kallistoOut = "abundance.tsv"){
    countMatrix <- readCountFile(tableCountPath,sepCharacter)
    designExperiment <- rep(groupName, each = numberReplics)
    result <- NULL

    result$bayseq<-runBaySeq(countMatrix,
                             designExperiment)

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

    # DESeq2 kallisto
    designExperimentDeseq2 <- colnames(countMatrix)
    if(typeof(countMatrix) == "double"){
        result$deseq2 <- runDeseq2(countMatrix,
                              groupName,
                              designExperimentDeseq2,
                              kallistoReport,
                              kallistoDir,
                              kallistoSubDir,
                              kallistoOut)
    }else{
        result$deseq2 <- runDeseq2(countMatrix,
                                   groupName,
                                   designExperimentDeseq2)
        # SAMSeq only count data
        print("**** SAMSeq run CANCELLED, enabled for count data only.")
        result$samseq<-runSamSeq(countMatrix,
                                 numberReplics,
                                 designExperiment)
    }
    return(result)
}
