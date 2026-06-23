#' Execute edgeR expression Analisys
#'
#' @param countMatrix either a matrix of raw (read) counts.
#' @param numberReplics number of replicate (technical or biologcal) integer
#' @param designExperiment replicate and treatment by samples
#' @param methodNorm normalization method to be used
#' @param methodAdjPvalue correction method, a character string. Can be abbreviated.
#' @param numberTopTable maximum number of genes to list
#'
#' @return limma report in data Frame fromat
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
#' toolResult <- NULL
#' toolResult$limma <- runLimma(counts, 2, rep(treats, each = 2))
runLimma <- function (countMatrix, numberReplics, designExperiment, methodNorm = "TMM", methodAdjPvalue = "BH", numberTopTable = 1000000){
    if (numberReplics <= 1){
        warning("limma-voom requires at least 2 replicates per condition. Skipping limma analysis.")
        return(NULL)
    }else {
        .check_package("edgeR", repo = "Bioconductor")
        .check_package("limma", repo = "Bioconductor")
        nf <- edgeR::calcNormFactors(countMatrix, method = methodNorm)
        condition <- factor(c(designExperiment))
        voom.data <- limma::voom(countMatrix, design = stats::model.matrix(~factor(condition)))
        voom.data$genes <- rownames(countMatrix)
        voom.fitlimma <- limma::lmFit(voom.data, design= stats::model.matrix(~factor(condition)))
        voom.fitbayes <- limma::eBayes(voom.fitlimma)
        voom.pvalues <- voom.fitbayes$p.value[, 2]
        voom.adjpvalues <- stats::p.adjust(voom.pvalues, method=methodAdjPvalue)
        data <- limma::topTable(voom.fitbayes, coef=ncol(designExperiment), number=numberTopTable)
        return(data)
    }
}
