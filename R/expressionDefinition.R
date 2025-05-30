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
#' @param pathOutput path to write output, need be a directory (default: ".")
#' #'
#' @export
#'
#' @examples
#' data(gse95077)
#' treats = c("BM", "JJ")
#' cons_result <- runExpression(numberReplics = 3,
#'                               groupName = treats,
#'                               rDataFrameCount = gse95077,
#'                               sepCharacter = ",",
#'                               experimentName = "test_cons",
#'                               controlDeseq2 = "BM",
#'                               contrastDeseq2 = "JJ",
#'                               outDirPath = "." )
#' expDef_result <- expressionDefinition(resultTool = cons_result, groups = treats)
expressionDefinition <- function(resultTool,groups = c(""),
                                 lfcMinLimma = -2,lfcMaxLimma = 2,pValueLimma = 0.05,FLimma = 0.8,
                                 lfcMinSamseq = -2,lfcMaxSamseq = 2,qValueSamseq = 0.05,scoreDSamseq = 0.8,
                                 lfcMinDeseq2 = -2,lfcMaxDeseq2 = 2,pValueDeseq2 = 0.05,
                                 lfcMinEdger = -2, lfcMaxEdger = 2, pValueEdger = 0.05,
                                 probNoiseq = 0.8,lfcMinKnowseq = -2,pValueKnowseq = 0.05,
                                 deClassEbseq = "DE", ppThresholdEbseq = 0.8,
                                 printResults = FALSE,pathOutput = "."){#, lfcMaxEbseq = 2, #lfcMinEbseq = -2){
    deList <- NULL
    if(!is.null(resultTool$edger)){ # edger
      deList$edger <- subset(resultTool$edger,(logFC <= lfcMinEdger  | `logFC` >= lfcMaxEdger) & `PValue` <= pValueEdger)
    }
    if(!is.null(resultTool$knowseq) && length(resultTool$knowseq) > 0){ # knowseq
      deList$knowseq <- subset(resultTool$knowseq,(`logFC` <= lfcMinKnowseq  | `logFC` >= lfcMaxKnowseq) & `P.Value` <= pValueKnowseq)
    }
    if(!is.null(resultTool$limma)){ #limma
      if(length(groups) > 2){
        deList$limma <- subset(resultTool$limma,(`F` >= FLimma) & (`P.Value` <= pValueLimma))
      }else{
        deList$limma <- subset(resultTool$limma,(`logFC` <= lfcMinLimma) | (`logFC` >= lfcMaxLimma) & (P.Value <= pValueLimma))
      }
    }
    if(!is.null(resultTool$noiseq)){ # NOISeq
      deList$noiseq <- subset(resultTool$noiseq,(prob >=probNoiseq))
    }
    if(!is.null(resultTool$ebseq)){ # ebseq
      ebseqDf <- as.data.frame(resultTool$ebseq)
      deList$ebseq <- subset(ebseqDf, resultTool$ebseq == deClassEbseq)
    }
    if(!is.null(resultTool$deseq2)){ # DESeq2
      deList$deseq2 <- subset(resultTool$deseq2,(`log2FoldChange` <= lfcMinDeseq2  | `log2FoldChange` >= lfcMaxDeseq2) & (pvalue <= pValueDeseq2))
    }
    if(!is.null(resultTool$samseq)){ #SAMSeq
      samseqDf <- as.data.frame(resultTool$samseq, row.names = NULL)
      if(length(groups) > 2){
        deList$samseq <- subset(samseqDf,(`Score(d)` >= scoreDSamseq) & (`q-value(%)` <= qValueSamseq))
      }else{
        deList$samseq <- subset(samseqDf, (`Fold Change` >= lfcMaxSamseq | `Fold Change` <= lfcMinLimma) & `q-value(%)` <= qValueSamseq)
      }
      row.names(deList$samseq) <- deList$samseq$`Gene ID`
    }
    if(printResults){
      tools <- names(deList)
      i <- 1
      for (deData in deList) {
        consexpressionR::writeResults(data = deData, toolName=paste0(tools[i],"_DE_"), pathOutput = pathOutput)
        i <- i + 1
      }
    }
    return(deList)
}
