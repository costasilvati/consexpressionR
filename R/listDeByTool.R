#' listDeByTool: creates a data.frame with genes (rows) and DE tools (columns).
#' Each tool that identifies a gene as DE is marked with 1.
#'
#' @param consexpressionList An \code{ExpressionResultSet} object with DE results.
#' @param geneNames Character vector. Names of all genes in the experiment.
#' @param deList List. Each element contains a character vector with DE gene names (from each tool).
#'
#' @return A data.frame with genes as rows and tools as columns; 1 indicates DE by that tool.
#' @export
#'
#' @examples
#' data(gse95077)
#' treats <- c("BM", "JJ")
#' cons_result <- runExpression(numberReplics = 3, groupName = treats, rDataFrameCount = gse95077, controlDeseq2 = "BM", contrastDeseq2 = "JJ" )
#' expDef_result <- expressionDefinition(resultTool = cons_result, groups = treats)
#' deByTool <- listDeByTool(cons_result, geneNames = rownames(gse95077), deList = expDef_result)
#'
#' #      DESeq2 limma edgeR
#' # gene1      1     1     1
#' # gene2      0     0     1
#' # gene3      1     0     1
#' # gene4      0     1     0
#' # gene5      0     0     0
#' # gene6      1     1     1
listDeByTool <- function(consexpressionList, geneNames, deList) {
  if (!inherits(consexpressionList, "ExpressionResultSet")) {
    stop("'consexpressionList' must be an ExpressionResultSet object.")
  }

  if (!is.character(geneNames)) {
    stop("'geneNames' must be a character vector.")
  }

  if (!is.list(deList)) {
    stop("'deList' must be a list.")
  }

  expressionList <- consexpressionList@results
  toolNames <- names(expressionList)

  deTool <- matrix(0L, nrow = length(geneNames), ncol = length(toolNames),
                   dimnames = list(geneNames, toolNames))

  for (i in seq_along(toolNames)) {
    degs <- rownames(deList[[i]])
    idx <- which(geneNames %in% degs)
    deTool[idx, i] <- 1L
  }

  return(as.data.frame(deTool))
}
