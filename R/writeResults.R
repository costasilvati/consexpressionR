#' Write a table or data.frame in defined output file
#'
#' @importFrom utils write.table
#' @param data table or data.frame dataset
#' @param toolName text with file name
#' @param sepCharacter character to separate data columns
#' @param pathOutput path to write output, needs to be a directory (default:".")
#'
#' @return The path of the file created (invisible)
#' @export
#'
#' @examples
#' df <- data.frame(col1 = c("treat1", "treat2", "treat3"),
#'                  col2 = c(1, 2, 3))
#' out <- writeResults(data = df, toolName = "test", pathOutput = tempdir())
#' file.exists(out)
writeResults <- function(data, toolName = "toolDE_x",
                         sepCharacter = "\t", pathOutput = ".") {
    if (!is.null(data)) {
        out <- createNameFileOutput(pathOutput, execName = toolName)
        utils::write.table(data, file = out, sep = sepCharacter, quote = FALSE)
        return(invisible(out))  # <- retorna o caminho sem poluir o console
    }
}
