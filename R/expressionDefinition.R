#' Identify DE genes by method-specific thresholds
#'
#' Applies specific filtering rules to each expression method result in the
#' ExpressionResultSet object and returns the DE genes per tool.
#'
#' @param resultTool An \code{ExpressionResultSet} object.
#' @param groups Character vector of group names or conditions.
#' @param lfcMinLimma,...,ppThresholdEbseq Numeric thresholds used per method (see details).
#' @param printResults Logical. Whether to write results to disk.
#' @param pathOutput Character. Directory path to save result tables.
#'
#' @return A list of data.frames with DE genes per method.
#' @export
#'
#' @examples
#' data(gse95077)
#' treats <- c("BM", "JJ")
#' res <- runExpression(numberReplics = 3, groupName = treats, rDataFrameCount = gse95077, controlDeseq2 = "BM", contrastDeseq2 = "JJ" )
#' expDef <- expressionDefinition(resultTool = res, groups = treats)
expressionDefinition <- function(resultTool, groups = c(""),
                                 lfcMinLimma = -2, lfcMaxLimma = 2, pValueLimma = 0.05, FLimma = 0.8,
                                 lfcMinSamseq = -2, lfcMaxSamseq = 2, qValueSamseq = 0.05, scoreDSamseq = 0.8,
                                 lfcMinDeseq2 = -2, lfcMaxDeseq2 = 2, pValueDeseq2 = 0.05,
                                 lfcMinEdger = -2, lfcMaxEdger = 2, pValueEdger = 0.05,
                                 probNoiseq = 0.8,
                                 lfcMaxKnowseq = 2, lfcMinKnowseq = -2, pValueKnowseq = 0.05,
                                 deClassEbseq = "DE", ppThresholdEbseq = 0.8,
                                 printResults = FALSE, pathOutput = ".") {

  if (!inherits(resultTool, "ExpressionResultSet")) {
    stop("'resultTool' must be an ExpressionResultSet object.")
  }

  deList <- list()

  if (!is.null(resultTool@results$edger)) {
    deList$edger <- subset(resultTool@results$edger,
                           (`logFC` <= lfcMinEdger | `logFC` >= lfcMaxEdger) &
                             (`PValue` <= pValueEdger))
  }

  if (!is.null(resultTool@results$knowseq) && length(resultTool@results$knowseq) > 0) {
    deList$knowseq <- subset(resultTool@results$knowseq,
                             (`logFC` <= lfcMinKnowseq | `logFC` >= lfcMaxKnowseq) &
                               (`P.Value` <= pValueKnowseq))
  }

  if (!is.null(resultTool@results$limma)) {
    if (length(groups) > 2) {
      deList$limma <- subset(resultTool@results$limma,
                             (`F` >= FLimma) & (`P.Value` <= pValueLimma))
    } else {
      deList$limma <- subset(resultTool@results$limma,
                             ((`logFC` <= lfcMinLimma) | (`logFC` >= lfcMaxLimma)) &
                               (`P.Value` <= pValueLimma))
    }
  }

  if (!is.null(resultTool@results$noiseq)) {
    deList$noiseq <- subset(resultTool@results$noiseq, prob >= probNoiseq)
  }

  if (!is.null(resultTool@results$ebseq)) {
    ebseqDf <- as.data.frame(resultTool@results$ebseq)
    deList$ebseq <- subset(ebseqDf, resultTool@results$ebseq == deClassEbseq)
  }

  if (!is.null(resultTool@results$deseq2)) {
    deList$deseq2 <- subset(resultTool@results$deseq2,
                            ((`log2FoldChange` <= lfcMinDeseq2) |
                               (`log2FoldChange` >= lfcMaxDeseq2)) &
                              (pvalue <= pValueDeseq2))
  }

  if (!is.null(resultTool@results$samseq)) {
    samseqDf <- as.data.frame(resultTool@results$samseq, row.names = NULL)
    if (length(groups) > 2) {
      deList$samseq <- subset(samseqDf,
                              (`Score(d)` >= scoreDSamseq) &
                                (`q-value(%)` <= qValueSamseq))
    } else {
      deList$samseq <- subset(samseqDf,
                              ((`Fold Change` >= lfcMaxSamseq) |
                                 (`Fold Change` <= lfcMinSamseq)) &
                                (`q-value(%)` <= qValueSamseq))
    }
    row.names(deList$samseq) <- deList$samseq$`Gene ID`
  }

  if (printResults && length(deList) > 0) {
    tools <- names(deList)
    for (i in seq_along(deList)) {
      consexpressionR::writeResults(
        data = deList[[i]],
        toolName = paste0(tools[i], "_DE_"),
        pathOutput = pathOutput
      )
    }
  }

  return(deList)
}
