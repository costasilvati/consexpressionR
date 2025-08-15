#' consensusList: make a list of genes considered DE by at least 'threshold' tools
#'
#' @param consexpressionList An \code{ExpressionResultSet} object containing the results from each analysis method.
#' @param deTool data.frame. Output of 'listDeByTool' function
#' @param threshold Integer. Minimum number of tools that must agree on a gene being DE (default: 5)
#'
#' @return An \code{ExpressionResultSet} object with the consensus slot populated
#' @export
#'
#' @examples
#' data(gse95077)
#' treats <- c("BM", "JJ")
#' cons_result <- runExpression(numberReplics = 3, groupName = treats, rDataFrameCount = gse95077, controlDeseq2 = "BM", contrastDeseq2 = "JJ" )
#' expDef_result <- expressionDefinition(resultTool = cons_result, groups = treats)
#' deByTool <- listDeByTool(cons_result, row.names(gse95077), expDef_result)
#' cons_result <- consensusList(cons_result, deByTool)
#'
#' summary(cons_result)
consensusList <- function(consexpressionList,
                          deTool,
                          threshold = 5) {
  if (!inherits(consexpressionList, "ExpressionResultSet")) {
    stop("'consexpressionList' must be an ExpressionResultSet object.")
  }

  if (!is.data.frame(deTool)) {
    stop("'deTool' must be a data.frame.")
  }

  deTool$nDE <- rowSums(deTool)
  consensus <- deTool$nDE >= threshold
  deCons <- subset(deTool, consensus)

  consensusList <- list()
  toolNames <- names(consexpressionList@results)

  for (i in seq_along(toolNames)) {
    tool <- toolNames[i]
    frame <- as.data.frame(consexpressionList@results[[tool]])
    itens <- frame[rownames(deCons), , drop = FALSE]
    consensusList[[tool]] <- itens
  }
  names(consensusList) <- toolNames

  consexpressionList@consensus <- consensusList
  return(consexpressionList)
}
