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
runEdger <- function (countMatrix,
                      numberReplics,
                      desingExperiment,
                      methNorm = "TMM"){
    group <- c(desingExperiment)
    y.dge <- edgeR::DGEList(counts = countMatrix, group = group)
    if (numberReplics < 1){
        print('Replicates not found by edgeR. EdgeR should be executed manual form.')
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
