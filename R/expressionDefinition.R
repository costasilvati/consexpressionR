#' Define who genes are deferentially expressed by method
#'
#' @param resultTool table results by consexpressionR function
#' @param groups text, name of samples or treatment
#' @param lfcMinLimma minimum value to consider of log Fold Change (default: -2).
#' @param lfcMaxLimma maximum value to consider of log Fold Change (default: 2)
#' @param pValueLimma minimum P-Value to consider (default: 0.05)
#' @param lfcMinSamseq minimum value to consider of log Fold Change (default: -2).
#' @param lfcMaxSamseq maximum value to consider of log Fold Change (default: 2)
#' @param qValueSamseq q-value is a measure of the statistical significance of the difference in expression between the compared groups, taking the problem of multiple comparisons into account (default: 0.8)
#' @param scoreDSamseq statistic used to evaluate the significance of each gene in multiclass analysis
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
# @param lfcMinEbseq minimum value to consider of log Fold Change (default: -2).
# @param lfcMaxEbseq maximum value to consider of log Fold Change (default: 2)
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
                                 qValueSamseq = 0.8,
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
                                 ppThresholdEbseq = 0.8){
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
      consexpressionR::writeResults(deList$limma,"limmaDE")
    }
    if(!is.null(resultTool$samseq)){ #SAMSeq
      samseqDf <- as.data.frame(resultTool$samseq, row.names = NULL)
      if(length(groups) > 2){
        deList$samseq <- dplyr::filter(samseqDf,
                                       ((`Score(d)` >= scoreDSamseq) & (`q-value(%)` <= qValueSamseq)))
      }else{
        deList$samseq <- dplyr::filter(samseqDf,
                                       ((`Fold Change` >= lfcMinSamseq  & `Fold Change` >= lfcMaxSamseq) & (`q-value(%)` <= qValueSamseq)))
      }

      row.names(deList$samseq) <- deList$samseq$`Gene ID`
      consexpressionR::writeResults(deList$samseq,"SAMSeqDE")
    }
    if(!is.null(resultTool$deseq2)){ # DESeq2
      deList$deseq2 <- dplyr::filter(resultTool$deseq2,
                                     ((log2FoldChange <= lfcMinDeseq2  | log2FoldChange >= lfcMaxDeseq2) & (pvalue <= pValueDeseq2)))
      consexpressionR::writeResults(deList$deseq2,"DESeq2DE")
    }
    if(!is.null(resultTool$edger)){ # edger
      deList$edger <- dplyr::filter(resultTool$edger,
                                    ((logFC <= lfcMinEdger  | logFC >= lfcMaxEdger) & `PValue` <= pValueEdger))
      consexpressionR::writeResults(deList$edger,"edgerDE")
    }
    if(!is.null(resultTool$noiseq)){ # NOISeq
      deList$noiseq <- dplyr::filter(resultTool$noiseq,
                                     (prob >=probNoiseq))
      consexpressionR::writeResults(deList$noiseq,
                                    "NOISeqDE")
    }
    if(!is.null(resultTool$knowseq)){ # knowseq
      deList$knowseq <- dplyr::filter(resultTool$knowseq,
                                      ((logFC <= lfcMinKnowseq  | logFC >= lfcMaxKnowseq) & `P.Value` <= pValueKnowseq))
      consexpressionR::writeResults(deList$knowseq,"knowseqDE")
    }
    if(!is.null(resultTool$ebseq)){ # ebseq
      ebseqDf <- as.data.frame(resultTool$ebseq)
      deList$ebseq <- dplyr::filter(ebseqDf,
                                      resultTool$ebseq == deClassEbseq)
      consexpressionR::writeResults(deList$ebseq,
                                      "EBSeqDE")
    }
    return(deList)
}
