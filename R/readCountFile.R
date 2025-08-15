#' Read table count File and creates a matrix from the given content file.
#' @importFrom utils read.csv
#' @param tableCountPath the name of the file which the data are to be read from. Each row of the table appears as one line of the file. If it does not contain an absolute path, the file name is relative to the current working directory, getwd(). This can be a compressed file (see file). Alternatively, file can be a readable text-mode connection (which will be opened for reading if necessary, and if so closed (and hence destroyed) at the end of the function call).
#' @param split the field separator character. Values on each line of the file are separated by this character. If sep = "" (the default for read.table) the separator is ‘,’, newlines or carriage returns.
#'
#' @return content of file in matrix R format
#' @export
#'
#' @examples
#' \donttest{
#' txt <- "gene,sample1,sample2\nGeneA,10,20\nGeneB,5,15"
#' con <- textConnection(txt)
#' mat <- readCountFile(tableCountPath = con, split = ",")
#' close(con)
#' }

readCountFile <- function(tableCountPath="data/table_count_df.csv",
                          split=","){
    tableCount <- utils::read.csv(tableCountPath,
                           sep=split,
                           row.names=1,
                           header=TRUE,
                           stringsAsFactors=FALSE,
                           na.strings = "NA")
    tableCount[is.na(tableCount)] <- 0
    return(tableCount)
}
