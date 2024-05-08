#' Define who genes are differentially expressed by method
#'
#' @param deTool table results by consexpressionR function
#' @param lfcMin minimum value to consider of log Fold Change (default: -2). Used by: KnowSeq, edgeR, limma, DESeq2 and SAMSeq.
#' @param lfcMax maximum value to consider of log Fold Change (default: 2), Used by: KnowSeq, edgeR, limma, DESeq2 and SAMSeq.
#' @param pValue maximum P-Value to consider (default: 0.05). Used by KnowSeq,, edgeR, limma and DESeq2, are consider p-Value >=.
#' @param probNoiseq floating point, minimum probability that a read count comes from a real gene expression peak, rather than background noiseq (default: 0.95)
#' @param qValue q-value is a measure of the statistical significance of the difference in expression between the compared groups, taking the problem of multiple comparisons into account (default: 1)
#' @param deClass name od class to consider by EBSeq (default: "DE")
#'
#' @return Object data.frame containing genes in the rows and tools executed in the columns, with a value of 1 for tools that considered the gene to be differentially expressed, and 0 for no DE.
#' @export
#'
#' @examples
#' cons_result <- consexpressionR(numberReplics = 3,
#'                               groupName = c("BM", "JJ"),
#'                               rDataFrameCount = table_count_df,
#'                               sepCharacter = ",",
#'                               experimentName = "test_cons",
#'                               outDirPath = "." )
#' expDef_result <- expressionDefinition(resultTool = cons_result)
expressionDefinition <- function(resultTool,
                                 lfcMin = -2,
                                 lfcMax = 2,
                                 pValue = 0.05,
                                 probNoiseq = 0.8,
                                 qValue = 0.8,
                                 deClass = "DE" ){
    deList <- NULL
    # if(!is.null(resultTool$bayseq)){ # baySeq
    #   deList$bayseq <- dplyr::filter(resultTool$bayseq, ((log2FoldChange <= lfcMin | log2FoldChange >= lfcMax) & (pvalue >=pValue)))
    #   consexpressionR::writeResults(deList$bayseq2De,"DESeq2DE")
    # }
    if(!is.null(resultTool$knowseq)){ # knowseq
      deList$knowseq <- dplyr::filter(resultTool$knowseq,
                                      ((logFC <= lfcMin  | logFC >= lfcMax) & `P.Value` <= pValue))
      consexpressionR::writeResults(deList$knowseq,"knowseqDE")
    }
    if(!is.null(resultTool$edger)){ # edger
        deList$edger <- dplyr::filter(resultTool$edger,
                                        ((logFC <= lfcMin  | logFC >= lfcMax) & `PValue` <= pValue))
        consexpressionR::writeResults(deList$edger,"edgerDE")
    }
    if(!is.null(resultTool$limma)){ #limma
        deList$limma <- dplyr::filter(resultTool$limma,
                                        ((logFC <= lfcMin | logFC >= lfcMax) & `P.Value` <= pValue))
        consexpressionR::writeResults(deList$limma,"limmaDE")
    }
    if(!is.null(resultTool$ebseq)){ # ebseq
        ebseqDf <- as.data.frame(resultTool$ebseq,
                                 row.names = NULL)
        deList$ebseq <- dplyr::filter(ebseqDf,
                                        resultTool$ebseq == deClass)
        consexpressionR::writeResults(deList$ebseq,
                                      "EBSeqDE")
    }
    if(!is.null(resultTool$noiseq)){ # NOISeq
        deList$noiseq <- dplyr::filter(resultTool$noiseq,
                                         (prob >=probNoiseq))
        consexpressionR::writeResults(deList$noiseq,
                                      "NOISeqDE")
    }
    if(!is.null(resultTool$deseq2)){ # DESeq2
        deList$deseq2 <- dplyr::filter(resultTool$deseq2,
                                         ((log2FoldChange <= lfcMin  | log2FoldChange >= lfcMax) & (pvalue <= pValue)))
        consexpressionR::writeResults(deList$deseq2,"DESeq2DE")
    }
    if(!is.null(resultTool$samseq)){ #SAMSeq
      samseqDf <- as.data.frame(resultTool$samseq, row.names = NULL)
      deList$samseq <- dplyr::filter(samseqDf,
                                       ((`Fold Change` >= lfcMin  & `Fold Change` >= lfcMax) & (`q-value(%)` <= qValue)))
      row.names(deList$samseq) <- deList$samseq$`Gene ID`
      consexpressionR::writeResults(deList$samseq,"SAMSeqDE")
    }
    return(deList)
}
