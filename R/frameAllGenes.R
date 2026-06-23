#' frameAllGenes
#' @description
#' This function creates a data.frame for each result of the diffrential expression method performed. It was created to
#' ensure that all genes been listed, even thought the gene does not appear as diffrentially expressed, in this case NA data is inserted.
#' This functis is executed into runExpression() function, so there is no need to call it separatelly.
#'
#' @param cons_list A named list of data.frames, each corresponding to a method (e.g., edgeR, DESeq2, etc.)
#' @param countMatrix RData object with table count, where line name is gene name, and name column is treat name
#'
#' @returns list with analysis expression containing: all methods executed and their results, ensuring thar all genes has been listed;
#' @export
#'
#' @examples
#' set.seed(42)
#' counts <- matrix(
#'   as.integer(c(
#'     rnbinom(200, mu = 10,  size = 1), 
#'     rnbinom(200, mu = 100, size = 1) 
#'   )),
#'   nrow = 100,
#'   dimnames = list(paste0("gene", seq_len(100)),
#'                     paste0("sample", seq_len(4)))
#' )
#' groups_info <- c("control", "control", "treated", "treated")
#' treats <- c("control", "treated")
#' cons_result <- runExpression(
#'   numberReplics = 2,
#'   groupName = treats,
#'   rDataFrameCount = counts,
#'   experimentName = "test_cons",
#'   controlDeseq2 = "control",
#'   contrastDeseq2 = "treated",
#'   outDirPath = tempdir(),
#'   printResults = FALSE
#' )
#'
#' genes_list <- results(cons_result)
#' framed <- frameAllGenes(cons_list = genes_list, countMatrix = counts)
#' head(framed)
frameAllGenes <- function(cons_list,countMatrix){
    newList <- list()
    tools <- names(cons_list)
    i <- 1
    for (tool_exp in cons_list) {
        result_frame <- as.data.frame(tool_exp)
        result_frame$gene_id <- rownames(result_frame)
        countMatrix_df <- as.data.frame(countMatrix)
        countMatrix_df$gene_id <- rownames(countMatrix_df)
        result <- merge(countMatrix_df,result_frame,by = "gene_id",all.x = TRUE)
        rownames(result) <- result$gene_id
        countMatrix_df$gene_id <- NULL
        result <- result[, !names(result) %in% names(countMatrix_df)]
        rownames(result) <- result$gene_id
        result$gene_id <- NULL
        newList[[i]] <- result
        message(" - ", i, tools[[i]])
        i <- i+1
    }
    names(newList) <- tools
    return(newList)
}
