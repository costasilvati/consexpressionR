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
