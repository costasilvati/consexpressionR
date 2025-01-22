frameAllGenes <- function(cons_list, countMatrix){
  newList <- list()
  tools <- names(cons_list)
  i <- 1
  for (tool_exp in cons_list) {
    result_frame <- as.data.frame(tool_exp)
    result_frame$gene_id <- rownames(result_frame)

    countMatrix_df <- as.data.frame(countMatrix)
    countMatrix_df$gene_id <- rownames(countMatrix_df)
    result <- merge(countMatrix_df,
                    result_frame,
                    by = "gene_id",
                    all.x = TRUE)
    rownames(result) <- result$gene_id
    countMatrix_df$gene_id <- NULL
    result <- result[, !names(result) %in% names(countMatrix_df)]
    rownames(result) <- result$gene_id
    result$gene_id <- NULL
    newList[[i]] <- result
    cat(" - ", i, tools[[i]])
    i <- i+1
  }
  names(newList) <- tools
  return(newList)
}
