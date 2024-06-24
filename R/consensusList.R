#' consensusList: make a list of genes was consider DE by threshold tools
#'
#' @param consexpressionList list. Output by 'consexpressionR' function
#' @param deTool data.frame. Output of 'deByTool' function
#' @param threshold integer. Number of tool was consider a gene as DE to filter (default: 5)
#'
#' @return data.frame with DE genes and values of limma analysis
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
#' consList <- consensusList(cons_result,deByTool)
consensusList <- function(consexpressionList,
                          deTool,
                          threshold = 5){
  deTool$nDE <- rowSums(deTool)
  consensus <- deTool$nDE >= threshold
  deCons <- subset(deTool, consensus)
  newList <- list()
  toolNames <- names(consexpressionList)
  for(i in 1:length(consexpressionList)){
    frame <- as.data.frame(consexpressionList[[i]])
    newList[[i]] <- frame[row.names(deCons),, drop= FALSE]
  }
  names(newList) <- toolNames
  return(newList)
}
