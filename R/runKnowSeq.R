runKnowSeq <- function(countMatrix,
                       groupName,
                       numberReplic,
                       filterId="ensembl_gene_id"){
  designExperiment <- rep(groupName, each = numberReplic)
  myAnnotation <- getGenesAnnotation(row.names(countMatrix),
                                     filter=filterId)
  designExperiment <- rep(groupName, each=numberReplic)
  knowSeq <- DEGsExtraction(countMatrix,
                            labels = designExperiment)
  return(knowSeq$DEG_Results$DEGs_Table)
}
