consensusList <- function(deTool,
                          threshold = 2){
  deTool$nDE <- rowSums(deTool)
  consensus <- deTool$nDE >= threshold
  deCons <- subset(deTool, consensus)
}
