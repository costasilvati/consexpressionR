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
#' set.seed(42)
#' counts <- matrix(
#'   as.integer(c(
#'     rnbinom(200, mu = 10, size = 1),
#'     rnbinom(200, mu = 100, size = 1)
#'   )),
#'   nrow = 100,
#'   dimnames = list(
#'     paste0("gene", seq_len(100)),
#'     paste0("sample", seq_len(4))
#'   )
#' )
#' groups_info <- c("control", "control", "treated", "treated")
#' treats <- c("control", "treated")
#' cons_result <- runExpression(2, treats,
#'   rDataFrameCount = counts,
#'   sepCharacter = ",", experimentName = "test_cons",
#'   outDirPath = "."
#' )
#' expDef_result <- expressionDefinition(resultTool = cons_result, groups = treats)
#' deByTool <- listDeByTool(
#'   consexpressionList = cons_result,
#'   geneNames = rownames(counts),
#'   deList = expDef_result
#' )
#' if (requireNamespace("UpSetR", quietly = TRUE)) {
#'   upSetPlotData <- upSetPlotTools(
#'     df = deByTool, condition = "Control_vs_Treat",
#'     pathOut = ".", writeData = FALSE
#'   )
#' }
upSetPlotTools <- function(df, condition = "condition_default", pathOut = ".", writeData = TRUE) {
  df[is.na(df)] <- 0
  .check_package("UpSetR", repo = "CRAN")
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    warning(
      "The 'UpSetR' package is not installed. ",
      "Install using R command: install.packages('UpSetR'). ",
      "This analysis will be skipped."
    )
    return(NULL)
  }
  meu_grafico <- UpSetR::upset(df, sets = colnames(df), sets.bar.color = "#56B4E9", order.by = "freq", empty.intersections = "on")
  df$number_of_methods <- rowSums(df)
  if (writeData) {
    grDevices::pdf(paste0(pathOut, "upset_plot_", condition, ".pdf"), width = 10, height = 7.5, onefile = FALSE)
    grDevices::dev.off()
    utils::write.csv(df, file = paste0(pathOut, "upset_plot_", condition, ".csv"))
  } else {
    UpSetR::upset(as.data.frame(df), sets = colnames(df), sets.bar.color = "#56B4E9", order.by = "freq", empty.intersections = "on")
  }
  return(df)
}
