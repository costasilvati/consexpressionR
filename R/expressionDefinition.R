#' Identify DE genes by method-specific thresholds
#'
#' Applies specific filtering rules to each expression method result in the
#' ExpressionResultSet object and returns the DE genes per tool.
#'
#' @param resultTool An \code{ExpressionResultSet} object.
#' @param groups Character vector of group names or conditions.
#' @param lfcMinLimma Minimum logFC threshold for limma.
#' @param lfcMaxLimma Maximum logFC threshold for limma.
#' @param pValueLimma P-value cutoff for limma.
#' @param FLimma Minimum F-statistic for multi-group limma.
#' @param lfcMinSamseq Minimum fold-change for SAMSeq.
#' @param lfcMaxSamseq Maximum fold-change for SAMSeq.
#' @param qValueSamseq q-value cutoff for SAMSeq.
#' @param scoreDSamseq Score(d) threshold for SAMSeq multi-class.
#' @param lfcMinDeseq2 Minimum log2 fold-change for DESeq2.
#' @param lfcMaxDeseq2 Maximum log2 fold-change for DESeq2.
#' @param pValueDeseq2 P-value cutoff for DESeq2.
#' @param lfcMinEdger Minimum logFC for edgeR.
#' @param lfcMaxEdger Maximum logFC for edgeR.
#' @param pValueEdger P-value cutoff for edgeR.
#' @param probNoiseq Probability threshold for NOISeq.
#' @param lfcMaxKnowseq Maximum logFC for KnowSeq.
#' @param lfcMinKnowseq Minimum logFC for KnowSeq.
#' @param pValueKnowseq P-value cutoff for KnowSeq.
#' @param deClassEbseq DE class to filter for EBSeq (e.g. "DE" or "EE").
#' @param ppThresholdEbseq Posterior probability threshold for EBSeq.
#' @param printResults Logical. Whether to write results to disk.
#' @param pathOutput Character. Directory path to save result tables.
#'
#' @return A list of data.frames with DE genes per method.
#' @export
#'
#' @examples
#' set.seed(42)
#' counts <- matrix(
#'   as.integer(c(
#'     rnbinom(200, mu = 10,  size = 1), 
#'     rnbinom(200, mu = 100, size = 1) 
#'   )),
#'   nrow = 100,
#'   dimnames = list(paste0("gene", seq_len(100)),
#'                     paste0("sample", seq_len(4)))
#' )
#' groups_info <- c("control", "control", "treated", "treated")
#' treats <- c("control", "treated")
#' res <- runExpression(numberReplics = 2,
#'                      groupName = treats,
#'                      rDataFrameCount = counts,
#'                      controlDeseq2 = "control",
#'                      contrastDeseq2 = "treated" )
#' expDef <- expressionDefinition(res, treats)
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
