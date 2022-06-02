runLimma <- function (countMatrix, numberReplics, designExperiment, limmaOutPath, methodNorm = "TMM", methodAdjPvalue = "BH", numberTopTable = 1000000){
    library("limma")
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