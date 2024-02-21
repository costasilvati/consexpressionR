#' Define who genes are differentially expressed by method
#'
#' @param deTool table results by consexpression2 function
#' @param lfc minimum value to consider of log Fold Change
#'
#' @return void
#' @export
#'
#' @examples
expressionDefinition <- function(resultTool,
                                 lfcMin = -2,
                                 lfcMax = 2,
                                 pValue = 0.05,
                                 prob = 0.95,
                                 qValue = 1,
                                 deClass = "DE" ){
    deList <- NULL
    # if(!is.null(resultTool$bayseq)){ # baySeq
    #   deList$bayseq <- dplyr::filter(resultTool$bayseq, ((log2FoldChange <= lfcMin | log2FoldChange >= lfcMax) & (pvalue >=pValue)))
    #   consexpression2::writeResults(deList$bayseq2De,"DESeq2DE")
    # }
    if(!is.null(resultTool$knowseq)){ # knowseq
      deList$knowseq <- dplyr::filter(resultTool$knowseq,
                                      ((logFC <= lfcMin  | logFC >= lfcMax) & `P.Value` <= pValue))
      consexpression2::writeResults(deList$knowseq,"knowseqDE")
    }
    if(!is.null(resultTool$edger)){ # edger
        deList$edger <- dplyr::filter(resultTool$edger,
                                        ((logFC <= lfcMin  | logFC >= lfcMax) & `PValue` <=pValue))
        consexpression2::writeResults(deList$edger,"edgerDE")
    }
    if(!is.null(resultTool$limma)){ #limma
        deList$limma <- dplyr::filter(resultTool$limma,
                                        ((logFC <= lfcMin | logFC >= lfcMax) & `P.Value` <= pValue))
        consexpression2::writeResults(deList$limma,"limmaDE")
    }
    if(!is.null(resultTool$ebseq)){ # ebseq
        ebseqDf <- as.data.frame(resultTool$ebseq,
                                 row.names = NULL)
        deList$ebseq <- dplyr::filter(ebseqDf,
                                        resultTool$ebseq == deClass)
        consexpression2::writeResults(deList$ebseq,
                                      "EBSeqDE")
    }
    if(!is.null(resultTool$noiseq)){ # NOISeq
        deList$noiseq <- dplyr::filter(resultTool$noiseq,
                                         (prob >=prob))
        consexpression2::writeResults(deList$noiseq,
                                      "NOISeqDE")
    }
    if(!is.null(resultTool$deseq2)){ # DESeq2
        deList$deseq2 <- dplyr::filter(resultTool$deseq2,
                                         ((log2FoldChange <= lfcMin  | log2FoldChange >= lfcMax) & (pvalue <= pValue)))
        consexpression2::writeResults(deList$deseq2,"DESeq2DE")
    }
    if(!is.null(resultTool$samseq)){ #SAMSeq
      samseqDf <- as.data.frame(resultTool$samseq, row.names = NULL)
      deList$samseq <- dplyr::filter(samseqDf,
                                       ((`Fold Change` >= lfcMin  & `Fold Change` >= lfcMax) & (`q-value(%)` <= qValue)))
      row.names(deList$samseq) <- deList$samseq$`Gene ID`
      consexpression2::writeResults(deList$samseq,"SAMSeqDE")
    }
    return(deList)
}
