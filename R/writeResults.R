#' Write a table or data.frame in defined output file
#'
#' @param data table or data.frame dataset
#' @param toolName text with file name
#' @param sepCharacter character to separe data columns
#'
#' @return void
#' @export
#'
#' @examples
#'
writeResults <- function(data,
                         toolName="toolDE_x",
                         sepCharacter="\t"){
  if(!is.null(data)){
    out <- createNameFileOutput(".",
                                execName=toolName)
    utils::write.table(data,
                       out,
                       sep=sepCharacter,
                       quote = FALSE)
  }
}
