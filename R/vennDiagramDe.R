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
  venn::venn(xDe, borders = FALSE, zcolor = "style", ellipse = TRUE)
  dev.off()
}
