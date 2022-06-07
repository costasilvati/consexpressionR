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
#' limmaResult<-runLimma(countMatrix, numberReplics, designExperiment, createNameFileOutput(outDirPath,experimentName,execName='limma'))
runLimma <- function (countMatrix, numberReplics, designExperiment, limmaOutPath, methodNorm = "TMM", methodAdjPvalue = "BH", numberTopTable = 1000000){
    if (numberReplics <= 1){
        print('ERROR: limma-voom require more than one replics.')
    }else {
        nf <- calcNormFactors(countMatrix, method = methodNorm)
        condition = factor(c(designExperiment))
        voom.data <- voom(countMatrix, design = model.matrix(~factor(condition)))
        voom.data$genes = rownames(countMatrix)
        voom.fitlimma = lmFit(voom.data, design=model.matrix(~factor(condition)))
        voom.fitbayes = eBayes(voom.fitlimma)
        voom.pvalues = voom.fitbayes$p.value[, 2]
        voom.adjpvalues = p.adjust(voom.pvalues, method=methodAdjPvalue)
        # design <- group
        data <- topTable(voom.fitbayes, coef=ncol(design), number=numberTopTable)
        write.table(data, file= limmaOutPath, sep = "\t", quote = FALSE)
        return(data)
    }
}
