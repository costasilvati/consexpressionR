#' UpSet Plot shows intersections by DE indications
#'
#' @param df data.frame: output by 'listDeByTool' function.
#' @param condition string: name of analysed experminet (default: "condition_default")
#' @param pathOut string: path to write a pdf file with upset Plot (default: ".")
#' @param writeData boolean: TRUE if want write a PDF file (default: TRUE)
#'
#' @return data.frame with row sum by input df
#' @export
#'
#' @examples
#' data(gse95077)
#' treats = c("BM", "JJ")
#' cons_result <- runExpression(numberReplics = 3, groupName = treats,rDataFrameCount = gse95077,
#'                             sepCharacter = ",",experimentName = "test_cons",outDirPath = "." )
#' expDef_result <- expressionDefinition(resultTool = cons_result, groups = treats)
#' deByTool <- listDeByTool(consexpressionList = cons_result, geneNames = row.names(gse95077), deList = expDef_result)
#' upSetPlotData <- upSetPlotTools(df = deByTool, condition = "Control_vs_Treat", pathOut = ".", writeData = FALSE)
upSetPlotTools <- function(df, condition = "condition_default", pathOut = ".", writeData = TRUE){
  df[is.na(df)] <- 0
  meu_grafico <- UpSetR::upset(df, sets = colnames(df), sets.bar.color = "#56B4E9",order.by = "freq", empty.intersections = "on")
  df$number_of_methods <- rowSums(df)
  if(writeData){
    grDevices::pdf(paste0(pathOut,"upset_plot_",condition,".pdf"), width = 10, height = 7.5, onefile = FALSE)
    #print(meu_grafico)
    grDevices::dev.off()
    utils::write.csv(df, file = paste0(pathOut,"upset_plot_", condition,".csv"))
  }else{
    UpSetR::upset(as.data.frame(df), sets = colnames(df), sets.bar.color = "#56B4E9", order.by = "freq", empty.intersections = "on")
  }
  return(df)
}

