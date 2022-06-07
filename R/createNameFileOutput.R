#' Build a file output name
#'
#' @param outDirPath the name of the file which the data are to be read from. Each row of the table appears as one line of the file. If it does not contain an absolute path, the file name is relative to the current working directory, getwd(). This can be a compressed file (see file). Alternatively, file can be a readable text-mode connection (which will be opened for reading if necessary, and if so closed (and hence destroyed) at the end of the function call).
#' @param experimentName the name of experiment analised by consexpression2
#' @param execName name of the tool run, can be a gene expression tool
#'
#' @return string name of path file output
#' @export
#'
#' @examples
#' nameFile <- createNameFileOutput ("../", "new_experiment", "test_exec")
createNameFileOutput <- function (outDirPath=".", experimentName="experiment", execName="exec"){
  if(str_ends(outDirPath,'/', negate = TRUE)){
    outDirPath <- paste(outDirPath, '/', sep = "")
  }
  outputFileName <- paste(outDirPath, experimentName, '_consexpression_', execName, '.csv', sep = "")
  return(outputFileName)
}
