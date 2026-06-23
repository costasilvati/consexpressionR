#' Execute edgeR expression Analisys
#'
#' @param countMatrix either a matrix of raw (read) counts.
#' @param numberReplics number of replicate (technical or biologcal) integer
#' @param desingExperiment replicate and treatment by samples
#' @param methNorm normalization method to be used in edgeR::calcNormFactors(), default: "TMM"
#'
#' @return edgeR report in data Frame
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
#' treats <- c("control", "treated")
#' numberReplicsModel <- 2
#' designExperimentModel <- rep(treats, each = numberReplicsModel)
#' toolResult <- NULL
#' toolResult$edger <- runEdger(countMatrix = counts, numberReplics = numberReplicsModel,
#'                                desingExperiment = designExperimentModel)
runEdger <- function (countMatrix, numberReplics,desingExperiment, methNorm = "TMM"){
    group <- c(desingExperiment)
    m <- as.matrix(countMatrix)
    .check_package("edgeR", repo = "Bioconductor")
    y.dge <- edgeR::DGEList(counts = m, group = group)
    if (numberReplics < 1){
        message('Replicates not found by edgeR. EdgeR wasn`t executed.')
    }else if(numberReplics == 1){
        bcv <- 0.2
        y.et <- edgeR::exactTest(y.dge, dispersion = bcv^2)
        y.tp <- edgeR::topTags(y.et, n = 100000)
        # utils::write.table(y.tp$table, edgerOutPath, sep = "\t", quote = FALSE)
    }else{
        y.dge <- edgeR::calcNormFactors(y.dge, method=methNorm)
        y.dge <- edgeR::estimateDisp(y.dge)
        y.dge <- edgeR::estimateCommonDisp(y.dge)
        y.et <- edgeR::exactTest(y.dge)
        y.tp <- edgeR::topTags(y.et, n = 100000)
        y.pvalues <- y.et$table$PValue
    }
    return(y.tp$table)
}
