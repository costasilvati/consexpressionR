#' Write a table or data.frame in defined output file
#'
#' @importFrom utils write.table
#' @param data table or data.frame dataset
#' @param toolName text with file name
#' @param sepCharacter character to separe data columns
#' @param pathOutput path to write output, need be a directory (default: ".")
#'
#' @return void
#' @export
#'
#' @examples
#' df <- data.frame(col1 = c("treat1", "treat2", "treat3"), col2 = c(1, 2, 3))
#' writeResults(data = df, toolName = "test")
writeResults <- function(data,
                         toolName="toolDE_x",
                         sepCharacter="\t",
                         pathOutput = "."){
  if(!is.null(data)){
    out <- createNameFileOutput(pathOutput,
                                execName=toolName)
    write.table(data,
                file = out,
                sep=sepCharacter,
                quote = FALSE)
  }
}
