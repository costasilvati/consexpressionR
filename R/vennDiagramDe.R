#' Make a venn diagram with DE genes by consexpression2 results
#'
#' @param deList Dataset eith consexpression2 results
#' @param outPath Path to write image
#'
#' @return void
#' @export
#'
#' @examples
vennDiagramDe <- function(deList, outPath){
  xDe <- NULL
  xDe$edgeR <- row.names(deList$edgerDe)
  xDe$limma <- row.names(deList$limmaDe)
  xDe$NOISeq <- row.names(deList$noiseqDe)
  xDe$ebSeq <- row.names(deList$ebseqDe)
  xDe$DESeq2 <- row.names(deList$deseq2De)
  if(is.null(deList$samseqDe)){
    xDe$SAMSeq <- row.names(deList$samseqDe)
  }
  pdf(file=outPath)
  venn::venn(xDe, borders = FALSE, zcolor = "style", ellipse = TRUE, ilabels = TRUE)
  dev.off()
}
