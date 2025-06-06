#' Make expression analysis of 7 tools, and return a list of results by each tool.
#'
#' @importFrom utils read.csv
#' @param numberReplics number of replicate (technical or biological) by sample
#' @param groupName text, name of samples or treatment
#' @param tableCountPath path to csv file that contains count data or abundance data (local)
#' @param sepCharacter character used to split csv data, can be comma or tab
#' @param rDataFrameCount RData object with table count, where line name is gene name, and name column is treat name
#' @param experimentName text, name of experiment
#' @param outDirPath path to write output, need be a directory (local)
#' @param methodNormLimma normalization method to be used (limma)
#' @param methodAdjPvalueLimma correction method, a character string. Can be abbreviated (limma)
#' @param numberTopTableLimma maximum number of genes to list (limma)
#' @param printResults logical variable: TRUE print report by each tool, FALSE print only consensus result
#' @param methodNormEdgeR normalization method to be used in edgeR::calcNormFactors(), default: "TMM"
#' @param kNoiseq Counts equal to 0 are replaced by k. By default, k = 0.5
#' @param factorNoiseq A string indicating the name of factor whose levels are the conditions to be compared.
#' @param lcNoiseq Length correction is done by dividing expression by length^lc. By default, lc = 0
#' @param replicatesNoiseq In this argument, the type of replicates to be used is defined: "technical", "biological" or "no" replicates. By default, "technical" replicates option is chosen.
#' @param normNoiseq Normalization method t can be one of "rpkm" (default), "uqua" (upper quartile), "tmm" (trimmed mean of M) or "n" (no normalization). Default: "rpkm".
#' @param condExpNoiseq A vector containing the two conditions to be compared by the differential expression algorithm (needed when the factor contains more than 2 different conditions).
#' @param respTypeSamseq Problem type: "Quantitative" for a continuous parameter; "Two class unpaired" for two classes with unpaired observations; "Survival" for censored survival outcome; "Multiclass": more than 2 groups; "Two class paired" for two classes with paired observations. Default: "Two class unpaired".
#' @param npermSamseq Number of permutations used to estimate false discovery rates. Default 100
#' @param fdrEbseq parameter used in EBTest function: fdr False Discovery Rate cutt off
#' @param maxRoundEbseq parameter used in EBTest function: Number of iterations. The default value is 50.
#' @param methodDeResultsEbseq parameter used in GetDEResults function: "robust" or "classic". Using the "robust" option, EBSeq is more robust to genes with outliers and genes with extremely small variances. Using the "classic" option, the results will be more comparable to those obtained by using the GetPPMat() function from earlier version (<= 1.7.0) of EBSeq
#' @param filterIdKnowseq The attribute used as filter to return the rest of the attributes.
#' @param notSapiensKnowseq A boolean value that indicates if the user wants the human annotation or another annotation available in BiomaRt. The possible not human dataset can be consulted by calling the following function: biomaRt::listDatasets(useMart("ensembl")).
#' @param fitTypeDeseq2 either "parametric", "local", "mean", or "glmGamPoi" for the type of fitting of dispersions to the mean intensity
#' @param controlDeseq2 group of samples that represents control in experiment, used by DESeq2; Default is "".
#' @param contrastDeseq2 specific reference level that you need compare (this name should be in groupName List)
#' @param deNovoAanalysis boolean value (TRUE if dataset don`t have a reference genome)
#' @param progressShiny shiny app element, used to show execution progress
#'
#' @return list with all analysis expression
#' @export
#'
#' @examples
#' data(gse95077)
#' treats = c("BM", "JJ")
#' cons_result <- runExpression(numberReplics = 3, groupName = treats, rDataFrameCount = gse95077,
#'                               sepCharacter = ",", experimentName = "test_cons", controlDeseq2 = "BM",
#'                               contrastDeseq2 = "JJ", outDirPath = "." )
#' expDef_result <- expressionDefinition(resultTool = cons_result, groups = treats)
runExpression <- function (numberReplics, groupName, tableCountPath = "data/gse95077.csv",
                             sepCharacter=",", rDataFrameCount = NULL, experimentName="genericExperiment",
                             outDirPath="../../consexpression2_results/", printResults=FALSE,
                             fitTypeDeseq2 = "local", controlDeseq2 = "", contrastDeseq2 = "",
                             methodNormLimma = "TMM", methodAdjPvalueLimma = "BH",
                             numberTopTableLimma = 1000000, filterIdKnowseq="ensembl_gene_id",
                             notSapiensKnowseq = FALSE, methodNormEdgeR = "TMM", normNoiseq = "rpkm",
                             kNoiseq = 0.5, factorNoiseq="Tissue", lcNoiseq = 0,
                             replicatesNoiseq = "technical", condExpNoiseq = c(""),
                             respTypeSamseq = "Two class unpaired",  npermSamseq = 100,
                             fdrEbseq=0.05, maxRoundEbseq = 50, methodDeResultsEbseq = "robust",
                             deNovoAanalysis = FALSE, progressShiny = NULL){
    if(!is.null(rDataFrameCount)){
      countMatrix <- as.matrix(rDataFrameCount)
    }else{ countMatrix <- as.matrix(readCountFile(tableCountPath,sepCharacter))}
    if(!is.null(progressShiny)){ progressShiny(detail = "Count matrix") }
    designExperiment <- rep(groupName, each = numberReplics)
    resultTool <- list()
    if(!is.null(progressShiny) && (!is.null(var))){ progressShiny(detail = "Executing edger...") }
    resultTool$edger <- runEdger(countMatrix, numberReplics, designExperiment, methNorm = methodNormEdgeR)
    if(!is.null(progressShiny)){ progressShiny(detail = "Executing KnowSeq...") }
    resultTool$knowseq <- NULL
    tryCatch({
      if(deNovoAanalysis){
        cat("\n ------------ KnowSeq is not running to deNovo analysis!\n")
      }else{
        resultTool$knowseq <- runKnowSeq(count = countMatrix, groupName = groupName, numberReplic = numberReplics,
                                         filterId = filterIdKnowseq, notSapiens = notSapiensKnowseq)
        if(!is.null(resultTool$knowseq)){
          cat("\n ------------ KnowSeq executed!\n")
        }else{ cat("\n ------------ KnowSeq multiclass not Implemented!\n") }
      }
    }, error = function(e) {
      message(paste("\n \n ===== ERROR: KnowSeq execution is failed === \n",e,"\n"))
    })
    if(!is.null(progressShiny)){ progressShiny(detail = "Executing limma...") }
    tryCatch({
      resultTool$limma<-runLimma(countMatrix, numberReplics, designExperiment,
                            methodNormLimma, methodAdjPvalueLimma, numberTopTableLimma)
      cat("\n ------------ limma executed!\n")
    }, error = function(e) {
      message(paste("\n \n ===== ERROR: limma execution is failed === \n",e,"\n"))
    })
    if(!is.null(progressShiny)){ progressShiny(detail = "Executing NOISeq...") }
    tryCatch({
      resultTool$noiseq<-runNoiSeq(countMatrix = countMatrix, groups = groupName, designExperiment= designExperiment,
                                   normParm = normNoiseq, kParam =kNoiseq, factorParam=factorNoiseq,
                                   lcParam = lcNoiseq, replicatesParam = replicatesNoiseq, condExp = condExpNoiseq)
      cat("\n ------------ NOISeq executed! \n")
    }, error = function(e) {
      message(paste("\n \n ===== ERROR: NOISeq execution is failed === \n",e,"\n"))
    })
    if(!is.null(progressShiny)){progressShiny(detail = "Executing EBSeq...")}
    tryCatch({
      resultTool$ebseq <- runEbseq(as.matrix(countMatrix), designExperiment, fdr = fdrEbseq,
                               maxRound = maxRoundEbseq, methodDeResults = methodDeResultsEbseq, groups = groupName)
      cat("\n ------------ EBSeq executed!\n")
    }, error = function(e) {
      message(paste("\n \n ===== ERROR: EBSeq execution is failed === \n",e,"\n"))
    })
    if(!is.null(progressShiny)){progressShiny(detail = "Executing DESeq2...")}
    if(typeof(countMatrix) == "double"){
        resultTool$deseq2 <- NULL
        cat("\n ------------ DESeq2 run CANCELLED, enabled for count data only.\n")
        cat("**** SAMSeq run CANCELLED, enabled for count data only.\n")
        if(!is.null(progressShiny)){ progressShiny(detail = "SAMSeq canceled...") }
    }else{
      tryCatch({
        resultTool$deseq2 <- runDeseq2(countMatrix, groupName, numberReplics, controlGroup = controlDeseq2,
                                       contrastGroup = contrastDeseq2, fitTypeParam = fitTypeDeseq2)
          cat("\n ------------ DESeq2 executed!\n")
      }, error = function(e) {
        message(paste("\n \n ===== ERROR: DESeq2 execution is failed === \n",e,"\n"))
      })
      if(!is.null(progressShiny)){ progressShiny(detail = "Executing SAMSeq...") }
      tryCatch({
        resultTool$samseq <- runSamSeq(countMatrix, designExperiment, respType = respTypeSamseq, numberPermutations = npermSamseq)
        cat("\n ------------ SAMSeq executed!\n")
      }, error = function(e) {
        message(paste("\n \n ===== ERROR: SAMSeq execution is failed === \n",e,"\n"))
      })
      if(!is.null(progressShiny)){ progressShiny(detail = "Writting results...")}
    }
    if(printResults){
      i<-1
      tryCatch({
        for (data in resultTool) {
          utils::write.csv(data, file = paste0(outDirPath, experimentName, "_", names(resultTool)[i], ".csv"))
          save(data, file = paste0(outDirPath, experimentName, "_", names(resultTool)[i], ".RData"))
          i <- i +1
        }
        resultTool <- frameAllGenes(resultTool, countMatrix)
        save(resultTool, file = paste0(outDirPath, experimentName, "_toolsResult.RData"))
      }, error = function(e) {
        message(paste("\n \n ===== ERROR: WRITE FILES execution is failed === \n",e,"\n"))
      })
    }
    if(!is.null(progressShiny)){ progressShiny(detail = "Complete!") }
    cat("\n ------------ DIFERENTIAL EXPRESSION ANALYSIS WAS COMPLETE!\n")
    return(resultTool)
}
