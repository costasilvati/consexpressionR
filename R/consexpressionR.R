#' Main function. Make expression analysis of 7 tools and print output by consenso
#'
#' @param numberReplics number of replicate (technical or biological) by sample
#' @param groupName text, name of samples or treatment
#' @param tableCountPath path to file csv that contains cout data or abundance data (local)
#' @param sepCharacter character used to split csv data, can be comma or tab
#' @param rDataFrameCount RData object with table count, where line name is gene name, and name column is treat name
#' @param experimentName text, name of experiment
#' @param outDirPath path to write output, need be a directory (local)
#' @param methodNorm normalization method to be used (limma)
#' @param methodAdjPvalue correction method, a character string. Can be abbreviated (limma)
#' @param numberTopTable maximum number of genes to list (limma)
#' @param printResults logical variable: TRUE print report by each tool, FALSE print only consensus result
#' @param kallistoReport path to kallisto report file _quant
#' @param kallistoDir name of kallistto results
#' @param kallistoSubDir name of directory in kallisto results
#' @param kallistoOut path to kallisto .tsv file
#'
#' @return list with all analisys expression
#' @export
#'
#' @examples
#' exp_result <- consexpressionR(numberReplics = 3, groupName = c("BM", "JJ"), tableCountPath = "data/GSE95077_filtred.csv", sepCharacter = ",", experimentName = "GSE95077", outDirPath = "." )
consexpressionR <- function (numberReplics,
                             groupName,
                             tableCountPath,
                             sepCharacter=",",
                             rDataFrameCount = NULL,
                             experimentName="genericExperiment",
                             outDirPath="consexpression2_results/",
                             methodNorm = "TMM",
                             methodAdjPvalue = "BH",
                             numberTopTable = 1000000,
                             printResults=FALSE,
                             kallistoReport = "report.txt",
                             kallistoDir = "kallisto_quant",
                             kallistoSubDir = "expermient_kallisto",
                             kallistoOut = "abundance.tsv",
                             baySeqRespType = "Multiclass",
                             knowSeqFilter = "ensembl_gene_id"){
    if(!is.null(rDataFrameCount)){
      countMatrix <- as.matrix(rDataFrameCount)
    }else{
      countMatrix <- as.matrix(readCountFile(tableCountPath,sepCharacter))
    }
    designExperiment <- rep(groupName, each = numberReplics)
    resultTool <- NULL
    # resultTool$bayseq<-runBaySeq(countMatrix,
    #                          groupName,
    #                          numberReplics)
    # cat("baySeq executed!\n")
    resultTool$knowseq <- runKnowSeq(as.matrix(countMatrix),
                                     groupName = groupName,
                                     numberReplic = numberReplics,
                                     filterId = knowSeqFilter)
    cat("knowSeq executed!\n")
    resultTool$edger<-runEdger(countMatrix,
                          numberReplics,
                          designExperiment)
    cat("edger executed!\n")
    resultTool$limma<-runLimma(countMatrix,
                          numberReplics,
                          designExperiment,
                          methodNorm,
                          methodAdjPvalue,
                          numberTopTable)
    cat("limma executed!\n")
    resultTool$noiseq<-runNoiSeq(countMatrix,
                                 designExperiment)
    cat("NOISeq executed!")
    resultTool$ebseq <- runEbseq(countMatrix,
                             designExperiment)
    cat("ebseq executed!\n")
    # DESeq2 kallisto
    if(typeof(countMatrix) == "double"){
        resultTool$deseq2 <- runDeseq2(countMatrix,
                              groupName,
                              numberReplics,
                              designExperiment,
                              kallistoReport,
                              kallistoDir,
                              kallistoSubDir,
                              kallistoOut)
        cat("DESeq2 executed!\n")
        # SAMSeq only count data
        cat("**** SAMSeq run CANCELLED, enabled for count data only.\n")
    }else{
        resultTool$deseq2 <- runDeseq2(countMatrix = countMatrix,
                                   groupName = groupName,
                                   numberReplics = numberReplics,
                                   designExperiment = designExperiment)
        cat("DESeq2 executed!\n")
        resultTool$samseq <- runSamSeq(countMatrix,
                                 designExperiment, respType = "Two class unpaired")
        cat("SAMSeq executed!\n")
    }
    if(printResults){
      i<-1
      for (data in resultTool) {
        write.csv(data, file = paste0(outDirPath,
                                      experimentName,
                                      "_",
                                      names(resultTool)[i],
                                      ".csv")
                  )
        save(data, file = paste0(outDirPath,
                                      experimentName,
                                      "_",
                                      names(resultTool)[i],
                                      ".RData")
        )
        i <- i +1
      }
      save(resultTool,
           file = paste0(outDirPath,
                         experimentName,
                         "_toolsResult.RData"))
    }
    return(resultTool)
}

#' SRA010153_consexpression2 <- consexpression2(numberReplics = 7, groupName = c("UHR", "Brain"), tableCountPath = "data/SRA010153_filtred.csv", sepCharacter = ",", experimentName = "SRA010153", outDirPath = "/Volumes/SD128/consexpression2_testesOutput/SRA010153/")
#' GSE95077_consexpression2 <- consexpression2(numberReplics = 3, groupName = c("BM", "JJ"), tableCountPath = "data/GSE95077_filtred.csv", sepCharacter = ",", experimentName = "GSE95077", outDirPath = "/Volumes/SD128/consexpression2_testesOutput/GSE95077/", printResults = TRUE)
