#' Define who genes are differentially expressed by method
#'
#' @param deTool table results by consexpression2 function
#' @param lfc minimum value to consider of log Fold Change
#'
#' @return void
#' @export
#'
#' @examples
expressionDefinition <- function(deTool, lfc){
    deList <- NULL
    if(!is.null(deTool$edger)){ # edger
        deList$edgerDe <- dplyr::filter(deTool$edger, ((logFC <= -2.0 | logFC >= 2.0) & (PValue >= 0.05)))
        consexpression2::writeResults(deList$edgerDe,"edgerDE")
    }
    if(!is.null(deTool$limma)){
        deList$limmaDe <- dplyr::filter(deTool$limma, ((logFC <= -2.0 | logFC >= 2.0) & (P.Value >= 0.05)))
        consexpression2::writeResults(deList$limmaDe,"limmaDE")
    }
    if(!is.null(deTool$ebseq)){ # ebseq
        ebseqDf <- as.data.frame(deTool$ebseq, row.names = NULL)
        deList$ebseqDe <- dplyr::filter(ebseqDf, deTool$ebseq == "DE")
        consexpression2::writeResults(deList$ebseqDe,"EBSeqDE")
    }
    if(!is.null(deTool$noiseq)){
        deList$noiseqDe <- dplyr::filter(deTool$noiseq, (prob >= 0.95))
        consexpression2::writeResults(deList$noiseqDe,"NOISeqDE")
    }
    if(!is.null(deTool$desq2)){
        resOrdered <- deTool$desq2[order(deTool$desq2$pvalue),]
        deseq2Df <- as.data.frame(resOrdered)
        deList$deseq2De <- dplyr::filter(deseq2Df, ((log2FoldChange <= -2.0 | log2FoldChange >= 2.0) & (pvalue >= 0.05)))
        consexpression2::writeResults(deList$deseq2De,"DESeq2DE")
    }

}
