runEdger <- function (countMatrix, numberReplics, desingExperiment, edgerOutPath){
    library("edgeR")
    group <- c(desingExperiment)
    y.dge <- DGEList(counts = countMatrix, group = group)
    if (numberReplics < 1){
        print('Replicates not found by edgeR. EdgeR should be executed manual form.')
    }else if(numberReplics == 1){
        bcv <- 0.2
        y.et <- exactTest(y.dge, dispersion = bcv^2)
        y.tp <- topTags(y.et, n = 100000)
        # y.pvalues <- y.et$table$PValue
        write.table(y.tp$table, edgerOutPath, sep = "\t", quote = FALSE)
    }else{
        y.dge <- calcNormFactors(y.dge)
        y.dge <- estimateDisp(y.dge)
        y.dge <- estimateCommonDisp(y.dge)
        y.et <- exactTest(y.dge)
        y.tp <- topTags(y.et, n = 100000)
        y.pvalues <- y.et$table$PValue
        write.table(y.tp$table, edgerOutPath, sep = "\t", quote = FALSE)
    }
    return(y.tp$table)
}