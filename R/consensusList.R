#' consensusList: make a list of genes was consider DE by threshold tools
#'
#' @param consexpressionList list. Output by 'consexpressionR' function
#' @param deTool data.frame. Output of 'deByTool' function
#' @param threshold integer. Number of tool was consider a gene as DE to filter (default: 5)
#'
#' @return data.frame with DE genes and values analysis
#' @export
#'
#' @examples
#' data(gse95077)
#' treats <- c("BM", "JJ")
#' cons_result <- runExpression(numberReplics = 3,groupName = treats,
#'                             rDataFrameCount = gse95077,
#'                             sepCharacter = ",",
#'                             experimentName = "test_cons",
#'                             controlDeseq2 = "BM",
#'                             contrastDeseq2 = "JJ",
#'                             outDirPath = ".")
#' expDef_result <- expressionDefinition(resultTool = cons_result,
#'                                       groups = treats)
#' deByTool <- listDeByTool(consexpressionList = cons_result,
#'                          geneNames = row.names(gse95077),
#'                          deList = expDef_result)
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
    itens <- frame[row.names(deCons),, drop= FALSE]
    newList[[i]] <- itens
  }
  names(newList) <- toolNames
  return(newList)
}
