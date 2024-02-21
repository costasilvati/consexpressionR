# SRA010153_DE_Tool <- listDeByTool(consexpressionList = SRA010153_consexpression2, deList = SRA010153_deList)
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
      } else {
        print("Termo não encontrado.")
      }
    }
    t <- t + 1
  }
  return(deTool)
}
