#' Execute edgeR expression Analisys
#'
#' @param countMatrix either a matrix of raw (read) counts.
#' @param numberReplics number of replicate (technical or biologcal) integer
#' @param designExperiment replicate and treatment by samples
#' @param limmaOutPath path to write limma report in csv file format
#' @param methodNorm normalization method to be used
#' @param methodAdjPvalue correction method, a character string. Can be abbreviated.
#' @param numberTopTable maximum number of genes to list
#'
#' @return limma report in data Frame fromat
#' @export
#'
#' @examples
#'
runLimma <- function (countMatrix,
                      numberReplics,
                      designExperiment,
                      methodNorm = "TMM",
                      methodAdjPvalue = "BH",
                      numberTopTable = 1000000){
    if (numberReplics <= 1){
        print('ERROR: limma-voom require more than one replics.')
    }else {
        nf <- edgeR::calcNormFactors(countMatrix, method = methodNorm)
        condition = factor(c(designExperiment))
        voom.data <- limma::voom(countMatrix, design = stats::model.matrix(~factor(condition)))
        voom.data$genes = rownames(countMatrix)
        voom.fitlimma = limma::lmFit(voom.data, design= stats::model.matrix(~factor(condition)))
        voom.fitbayes = limma::eBayes(voom.fitlimma)
        voom.pvalues = voom.fitbayes$p.value[, 2]
        voom.adjpvalues = stats::p.adjust(voom.pvalues, method=methodAdjPvalue)
        # design <- group
        data <- limma::topTable(voom.fitbayes, coef=ncol(designExperiment), number=numberTopTable)
        return(data)
    }
}
