#' Makes a heatmap plot with differential expressed genes
#'
#' @param countData
#' @param designExpreiment
#' @param pdfOutName string: path + name file, to write a pdf file with Plot (default: "heatMap.pdf")
#' @param condition string: name of analysed experminet (default: "condition_default")
#' @param writePdf logical variable: TRUE print report by each tool, FALSE print only consensus result
#'
#' @return
#' @export
#'
#' @examples
heatMapPlot <- function(countData,
                        designExpreiment,
                        pdfOutName = "heatMap.pdf",
                        condition = "default",
                        writePdf = FALSE){

  col <- gplots::bluered(nrow(countData))
  sideColors <- gplots::bluered(ncol(countData))
  if(writePdf){
    grDevices::pdf(pdfOutName,
                   width = 15,
                   height = 12,
                   onefile = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  par(mar = c(5, 4, 4, 2) + 0.1)
  colnames(countData) <- designExpreiment
  gplots::heatmap.2(log2(as.matrix(countData)),
                    col = col,
                    scale = "row",
                    ColSideColors = sideColors,
                    key = TRUE,
                    symkey = FALSE,
                    density.info = "none",
                    cexRow = 0.8,
                    cexCol = 0.9,
                    keysize = 1.5,
                    margins = c(5, 7.5),
                    trace = "none",
                    srtCol = 25)
}
