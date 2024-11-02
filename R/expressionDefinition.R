#' Define who genes are deferentially expressed by method
#'
#' @param resultTool table results by consexpressionR function
#' @param groups text, name of samples or treatment
#' @param lfcMinLimma minimum value to consider of log Fold Change (default: -2).
#' @param lfcMaxLimma maximum value to consider of log Fold Change (default: 2)
#' @param pValueLimma minimum P-Value to consider (default: 0.05)
#' @param FLimma The F statistics is used to test the null hypotheses of all groups have the same median of expression
#' @param lfcMinSamseq minimum value to consider of log Fold Change (default: -2).
#' @param lfcMaxSamseq maximum value to consider of log Fold Change (default: 2)
#' @param qValueSamseq q-value is a measure of the statistical significance of the difference in expression between the compared groups, taking the problem of multiple comparisons into account (default: 0.8)
#' @param scoreDSamseq statistic used to evaluate the significance of each gene in multi class analysis
#' @param lfcMinDeseq2 minimum value to consider of log Fold Change (default: -2).
#' @param lfcMaxDeseq2 maximum value to consider of log Fold Change (default: 2)
#' @param pValueDeseq2 minimum P-Value to consider (default: 0.05)
#' @param lfcMinEdger minimum value to consider of log Fold Change (default: -2).
#' @param lfcMaxEdger maximum value to consider of log Fold Change (default: 2)
#' @param pValueEdger minimum P-Value to consider (default: 0.05)
#' @param probNoiseq floating point, minimum probability that a read count comes from a real gene expression peak, rather than background noiseq (default: 0.95)
#' @param lfcMinKnowseq minimum value to consider of log Fold Change (default: -2).
#' @param lfcMaxKnowseq maximum value to consider of log Fold Change (default: 2)
#' @param pValueKnowseq minimum P-Value to consider (default: 0.05)
#' @param deClassEbseq name of class to consider by EBSeq (default: "DE")
#' @param ppThresholdEbseq description
#' @param printResults logical variable: TRUE print report by each tool, FALSE print only consensus result
#' @export
#' @examples
#' library(cqn)
#' cons_result <- runExpression(numberReplics = 3,
#'                               groupName = c("BM", "JJ"),
#'                               rDataFrameCount = gse95077,
#'                               sepCharacter = ",",
#'                               experimentName = "test_cons",
#'                               outDirPath = "." )
#' expDef_result <- expressionDefinition(resultTool = cons_result)
expressionDefinition <- function(resultTool,
                                 groups = c(""),
                                 lfcMinLimma = -2,
                                 lfcMaxLimma = 2,
                                 pValueLimma = 0.05,
                                 FLimma = 0.8,
                                 lfcMinSamseq = -2,
                                 lfcMaxSamseq = 2,
                                 qValueSamseq = 0.05,
                                 lfcMinDeseq2 = -2,
                                 scoreDSamseq = 0.8,
                                 lfcMaxDeseq2 = 2,
                                 pValueDeseq2 = 0.05,
                                 lfcMinEdger = -2,
                                 lfcMaxEdger = 2,
                                 pValueEdger = 0.05,
                                 probNoiseq = 0.8,
                                 lfcMinKnowseq = -2,
                                 lfcMaxKnowseq = 2,
                                 pValueKnowseq = 0.05,
                                 deClassEbseq = "DE",
                                 ppThresholdEbseq = 0.8,
                                 printResults = FALSE){
                                 #, lfcMaxEbseq = 2,
                                 #lfcMinEbseq = -2){
    deList <- NULL
    if(!is.null(resultTool$limma)){ #limma
      if(length(groups) > 2){
        deList$limma <- dplyr::filter(resultTool$limma,
                                      ((F >= FLimma) & `P.Value` <= pValueLimma))
      }else{
        deList$limma <- dplyr::filter(resultTool$limma,
                                      ((logFC <= lfcMinLimma | logFC >= lfcMaxLimma) & `P.Value` <= pValueLimma))
      }
    }
    if(!is.null(resultTool$samseq)){ #SAMSeq
      samseqDf <- as.data.frame(resultTool$samseq, row.names = NULL)
      if(length(groups) > 2){
        deList$samseq <- dplyr::filter(samseqDf,
                                       ((`Score(d)` >= scoreDSamseq) & (`q-value(%)` <= qValueSamseq)))
      }else{
        deList$samseq <- dplyr::filter(samseqDf,
                                       ((`Fold Change` >= 2.0) | (`Fold Change` <= -2.0)) & `q-value(%)` <= 0.05)
      }
      row.names(deList$samseq) <- deList$samseq$`Gene ID`
    }
    if(!is.null(resultTool$deseq2)){ # DESeq2
      deList$deseq2 <- dplyr::filter(resultTool$deseq2,
                                     ((log2FoldChange <= lfcMinDeseq2  | log2FoldChange >= lfcMaxDeseq2) & (pvalue <= pValueDeseq2)))
    }
    if(!is.null(resultTool$edger)){ # edger
      deList$edger <- dplyr::filter(resultTool$edger,
                                    ((logFC <= lfcMinEdger  | logFC >= lfcMaxEdger) & `PValue` <= pValueEdger))
    }
    if(!is.null(resultTool$noiseq)){ # NOISeq
      deList$noiseq <- dplyr::filter(resultTool$noiseq,
                                     (prob >=probNoiseq))
    }
    if(!is.null(resultTool$knowseq)){ # knowseq
      deList$knowseq <- dplyr::filter(resultTool$knowseq,
                                      ((logFC <= lfcMinKnowseq  | logFC >= lfcMaxKnowseq) & `P.Value` <= pValueKnowseq))
    }
    if(!is.null(resultTool$ebseq)){ # ebseq
      ebseqDf <- as.data.frame(resultTool$ebseq)
      deList$ebseq <- dplyr::filter(ebseqDf,
                                      resultTool$ebseq == deClassEbseq)
    }

    if(printResults){
      tools <- names(deList)
      i <- 1
      for (deData in deList) {
        consexpressionR::writeResults(deData,paste0(tools[i], "DE"))
        i <- i + 1
      }
    }
    return(deList)
}
