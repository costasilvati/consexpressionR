#' listDeByTool: makes a data.frame with genes (lines) and all executed DE tool.
#' Each tool that identify gene as DE, are marked by 1.
#'
#' @param consexpressionList A list generated by consexpressionR function
#' @param deList List with names of all genes analysed in experiment
#'
#' @return data.frame with genes (lines) and DE tools (column), where 1 shows a
#' DE indication by tool.
#' @export
#'
#' @examples
#' library(cqn)
#' cons_result <- runExpression(numberReplics = 3, groupName = c("BM", "JJ"),
#'                               rDataFrameCount = gse95077,
#'                               sepCharacter = ",",
#'                               experimentName = "test_cons",
#'                               outDirPath = "." )
#' expDef_result <- expressionDefinition(resultTool = cons_result)
#' deByTool <- listDeByTool(cons_result, expDef_result)
listDeByTool <- function(consexpressionList, deList){
  deTool <- data.frame(matrix(0, ncol = length(consexpressionList),
                              nrow = nrow(consexpressionList[[2]])))
  colnames(deTool) <- names(consexpressionList)
  row.names(deTool) <- row.names(consexpressionList[[2]])
  toolNames <- names(consexpressionList)
  geneNames <- row.names(consexpressionList[[2]])
  t <- 1
  for (toolDe in consexpressionList) {
    genesDe <- rownames(toolDe)
    for (i in 1:length(geneNames)) {
      if(geneNames[i] %in% row.names(deList[[t]])){
        deTool[geneNames[i], t] <- 1
      }
    }
    t <- t + 1
  }
  return(deTool)
}
