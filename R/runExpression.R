#' Make expression analysis of 7 tools, and return a list of results by each tool.
#'
#' @importFrom utils read.csv
#' @param numberReplics number of replicate (technical or biological) by sample
#' @param groupName text, name of samples or treatment
#' @param tableCountPath path to file csv that contains cout data or abundance data (local)
#' @param sepCharacter character used to split csv data, can be comma or tab
#' @param rDataFrameCount RData object with table count, where line name is gene name, and name column is treat name
#' @param experimentName text, name of experiment
#' @param outDirPath path to write output, need be a directory (local)
#' @param methodNormLimma normalization method to be used (limma)
#' @param methodAdjPvalueLimma correction method, a character string. Can be abbreviated (limma)
#' @param numberTopTableLimma maximum number of genes to list (limma)
#' @param printResults logical variable: TRUE print report by each tool, FALSE print only consensus result
#' @param kallistoReport path to kallisto report file _quant
#' @param kallistoDir name of kallistto results
#' @param kallistoSubDir name of directory in kallisto results
#' @param kallistoOut path to kallisto .tsv file
#' @param methodNormEdgeR normalization method to be used in edgeR::calcNormFactors(), default: "TMM"
#' @param kNoiseq Counts equal to 0 are replaced by k. By default, k = 0.5
#' @param factorNoiseq A string indicating the name of factor whose levels are the conditions to be compared.
#' @param lcNoiseq Length correction is done by dividing expression by length^lc. By default, lc = 0
#' @param replicatesNoiseq In this argument, the type of replicates to be used is defined: "technical", "biological" or "no" replicates. By default, "technical" replicates option is chosen.
#' @param normNoiseq Normalization method t can be one of "rpkm" (default), "uqua" (upper quartile), "tmm" (trimmed mean of M) or "n" (no normalization). Default: "rpkm".
#' @param respTypeSamseq Problem type: "Quantitative" for a continuous parameter; "Two class unpaired" for two classes with unpaired observations; "Survival" for censored survival outcome; "Multiclass": more than 2 groups; "Two class paired" for two classes with paired observations. Default: "Two class unpaired".
#' @param npermSamseq Number of permutations used to estimate false discovery rates. Default 100
#' @param fdrEbseq parameter used in EBTest function: fdr False Discovery Rate cutt off
#' @param maxRoundEbseq parameter used in EBTest function: Number of iterations. The default value is 50.
#' @param methodDeResultsEbseq parameter used in GetDEResults function: "robust" or "classic". Using the "robust" option, EBSeq is more robust to genes with outliers and genes with extremely small variances. Using the "classic" option, the results will be more comparable to those obtained by using the GetPPMat() function from earlier version (<= 1.7.0) of EBSeq
#' @param filterIdKnowseq The attribute used as filter to return the rest of the attributes.
#' @param notSapiensKnowseq A boolean value that indicates if the user wants the human annotation or another annotation available in BiomaRt. The possible not human dataset can be consulted by calling the following function: biomaRt::listDatasets(useMart("ensembl")).
#' @param fitTypeDeseq2 either "parametric", "local", "mean", or "glmGamPoi" for the type of fitting of dispersions to the mean intensity
#' @return list with all analisys expression
#' @export
#'
#' @examples
#' cons_result <- runExpression(numberReplics = 3, groupName = c("BM", "JJ"),
#'                               rDataFrameCount = table_count_df,
#'                               sepCharacter = ",",
#'                               experimentName = "test_cons",
#'                               outDirPath = "." )
runExpression <- function (numberReplics,
                             groupName,
                             tableCountPath,
                             sepCharacter=",",
                             rDataFrameCount = NULL,
                             experimentName="genericExperiment",
                             outDirPath="consexpression2_results/",
                             printResults=FALSE,
                             fitTypeDeseq2 = "local",
                             kallistoReport = "report.txt",
                             kallistoDir = "kallisto_quant",
                             kallistoSubDir = "expermient_kallisto",
                             kallistoOut = "abundance.tsv",
                             methodNormLimma = "TMM",
                             methodAdjPvalueLimma = "BH",
                             numberTopTableLimma = 1000000,
                             filterIdKnowseq="ensembl_gene_id",
                             notSapiensKnowseq = FALSE,
                             methodNormEdgeR = "TMM",
                             normNoiseq = "rpkm",
                             kNoiseq = 0.5,
                             factorNoiseq="Tissue",
                             lcNoiseq = 0,
                             replicatesNoiseq = "technical",
                             respTypeSamseq = "Two class unpaired",
                             npermSamseq = 100,
                             fdrEbseq=0.05,
                             maxRoundEbseq = 50,
                             methodDeResultsEbseq = "robust"
                             ){
    if(!is.null(rDataFrameCount)){
      countMatrix <- as.matrix(rDataFrameCount)
    }else{
      countMatrix <- as.matrix(readCountFile(tableCountPath,sepCharacter))
    }
    designExperiment <- rep(groupName, each = numberReplics)
    resultTool <- NULL
    resultTool$knowseq <- runKnowSeq(as.matrix(countMatrix),
                                     groupName = groupName,
                                     numberReplic = numberReplics,
                                     filterId = filterIdKnowseq,
                                     notSapiens = notSapiensKnowseq)
    cat("knowSeq executed!\n")
    resultTool$edger<-runEdger(countMatrix,
                          numberReplics,
                          designExperiment,
                          methNorm = methodNormEdgeR)
    cat("edger executed!\n")
    resultTool$limma<-runLimma(countMatrix,
                          numberReplics,
                          designExperiment,
                          methodNormLimma,
                          methodAdjPvalueLimma,
                          numberTopTableLimma)
    cat("limma executed!\n")
    resultTool$noiseq<-runNoiSeq(countMatrix = countMatrix,
                                 designExperiment= designExperiment,
                                 normParm = normNoiseq,
                                 kParam =kNoiseq,
                                 factorParam=factorNoiseq,
                                 lcParam = lcNoiseq,
                                 replicatesParam = replicatesNoiseq)
    cat("NOISeq executed!")
    resultTool$ebseq <- runEbseq(countMatrix,
                             designExperiment,
                             fdr = fdrEbseq,
                             maxRound = maxRoundEbseq,
                             methodDeResults = methodDeResultsEbseq)
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
                              kallistoOut,
                              fitTypeParm = fitTypeDeseq2)
        cat("DESeq2 executed!\n")
        # SAMSeq only count data
        cat("**** SAMSeq run CANCELLED, enabled for count data only.\n")
    }else{
      resultTool$deseq2 <- runDeseq2(countMatrix,
                                     groupName,
                                     numberReplics,
                                     designExperiment,
                                     fitTypeParm = fitTypeDeseq2)
        cat("DESeq2 executed!\n")
        resultTool$samseq <- runSamSeq(countMatrix,
                                 designExperiment,
                                 respType = respTypeSamseq,
                                numberPermutations = npermSamseq)
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
    cat("DIFERENTIAL EXPRESSION ANALYSIS WAS COMPLETE!")
    return(resultTool)
}

#' SRA010153_consexpression2 <- consexpression2(numberReplics = 7, groupName = c("UHR", "Brain"), tableCountPath = "data/SRA010153_filtred.csv", sepCharacter = ",", experimentName = "SRA010153", outDirPath = "/Volumes/SD128/consexpression2_testesOutput/SRA010153/")
#' GSE95077_consexpression2 <- consexpression2(numberReplics = 3, groupName = c("BM", "JJ"), tableCountPath = "data/GSE95077_filtred.csv", sepCharacter = ",", experimentName = "GSE95077", outDirPath = "/Volumes/SD128/consexpression2_testesOutput/GSE95077/", printResults = TRUE)
