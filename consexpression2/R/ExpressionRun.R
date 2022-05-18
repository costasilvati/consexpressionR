# Title     : Run baySeq
# Objective : Execute baySeq with default paraemters by a count table
# Created by: julianacostasilva
# Created on: 30/03/22

# install.packages("styler") # formata
# https://cran.r-project.org/doc/manuals/R-exts.html#Creating-R-packages
# http://bioconductor.org/developers/package-submission/
# https://www.bioconductor.org/developers/how-to/coding-style/


readCountFile <- function (tableCountPath, split){
    tableCount <- read.csv(tableCountPath, sep=split, row.names=1, header=TRUE, stringsAsFactors=FALSE)
    return(as.matrix(tableCount))
}

createNameFileOutput <- function (outDirPath, experimentName, execName){
    library(stringr)
    if(str_ends(outDirPath,'/', negate = TRUE)){
        outDirPath <- paste(outDirPath, '/', sep = "")
    }
    outputFileName <- paste(outDirPath, experimentName, '_consexpression_', execName, '.csv', sep = "")
    return(outputFileName)
}

runBaySeq <- function  (countMatrix, designExperiment, out){
    library("baySeq")
    if(require("parallel")) cl <- makeCluster(4) else cl <- NULL
    groups <- list(NDE = designExperiment, DE = designExperiment)
    CD <- new("countData", data = countMatrix, replicates = designExperiment, groups = groups)
    libsizes(CD) <- getLibsizes(CD)
    CD <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl, equalDispersions = TRUE)
    CD <- getLikelihoods(CD, prs=c(0.5, 0.5), pET="BIC", cl=cl)
    write.table(topCounts(CD, group = "DE", number = 65000, normaliseData = TRUE), out, sep="\t", quote = FALSE)
    return (CD)
}

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

runNOISeq <- function (countMatrix, designExperiment, noiseqOutPtah){
    library("NOISeq")
    myfactors = data.frame(Tissue=c(designExperiment))
    mydata <- readData(data=countMatrix, factors=myfactors)
    # k?? lc?? factor??
    mynoiseq = noiseq(mydata, k=0.5, factor="Tissue", lc=1, replicates="technical")
    result <- mynoiseq@results[[1]]
    write.csv(result, file=noiseqOutPtah, sep="\t", quote=FALSE)
    return(mynoiseq)
}

installDependencies <- function (){ #SAMSeq install
    install.packages(c("matrixStats", "GSA", "shiny", "openxlsx", "Rcpp"))
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("impute")
    install.packages("devtools") if devtools not installed
    #library(devtools)
    #install_github("cran/samr")
    install.packages("samr")
}


runSamSeq <- function (countMatrix, split, designExperiment, replic, samseqOutPtah){
    library(samr)
    SAMseq.test <- SAMseq(countMatrix,
                          as.factor(designExperiment),
                          each=replic,
                          resp.type='Multiclass',
                          geneid=row.names(m),
                          genenames=row.names(m),
                          nperms=100)
    SAMseq.result.table <- rbind(SAMseq.test$siggenes.table$genes.up, SAMseq.test$siggenes.table$genes.lo)
    SAMseq.score <- rep(0, nrow(countMatrix))
    SAMseq.score[match(SAMseq.result.table[,1], rownames(countMatrix))]=as.numeric(SAMseq.result.table[,3])
    SAMseq.FDR = rep(1, nrow(countMatrix))
    SAMseq.FDR[match(SAMseq.result.table[,1], rownames(countMatrix))] = as.numeric(SAMseq.result.table[,5])/100
    write.table(SAMseq.result.table, file=samseqOutPtah, sep=split, quote=FALSE)
}

consexpression <- function (numberReplics=0, groupName,tableCountPath, split=",", experimentName='genericExperiment', outDirPath="."){
    # validações?
    countMatrix <- readCountFile(tableCountPath, split)
    designExperiment <- rep(groupName, each = numberReplics)
    # baySeq ok
    bayseqResult <- runBaySeq(countMatrix, replicates, createNameFileOutput(outDirPath,experimentName,execName = 'baySeq'))
    # edgeR ok
    edgerResult <- runEdger(countMatrix, numberReplics, designExperiment, createNameFileOutput(outDirPath,experimentName,execName = 'edgeR'))
    # limma ok
    limmaResult <- runLimma(countMatrix, numberReplics, designExperiment, createNameFileOutput(outDirPath,experimentName,execName = 'limma'))
    # NOISeq ok
    # noiseqResult <- runNOISeq(countMatrix, designExperiment, createNameFileOutput(outDirPath,experimentName,execName = 'NOISeq'))
}


numberReplics <- 3
groupName <- c("0b", "1b")
tableCountPath<- '/Users/julianacostasilva/OneDrive/ProjetoDoutorado/bioconvergencia/reads_RNApa/kallisto_quant_align_apa_1B_0B/gonoda/gonoda_tpm.csv'
split <- "\t"
experimentName <- 'RNApa_apa_1B_0B_gonoda'
outDirPath <- '/Volumes/SD128/bioconvergencia/reads_RNApa/'

#consexpression(numberReplics, groupName, tableCountPath,split,experimentName, outDirPath)

# SAMSeq ??
# samseqResult <- runSamseq(countMatrix, numberReplics, designExperiment, createNameFileOutput(outDirPath,experimentName,execName = 'SAMSeq'))
# EBSeq ??s
# ebseqResult <- runEbseq(countMatrix, numberReplics, designExperiment, createNameFileOutput(outDirPath,experimentName,execName = 'EBSeq'))
# DESeq ??
# deseqResult <- runDeseq(countMatrix, numberReplics, designExperiment, createNameFileOutput(outDirPath,experimentName,execName = 'DESeq'))
# DESeq2 ??
# deseq2OutPath = '/Volumes/SD128/bioconvergencia/reads_RNApa/RNApa_apa_1B_0B_gonoda_consexpression_deseq2.csv'
# deseq2Result <- runDeseq(countMatrix, numberReplics, designExperiment, deseq2OutPath)

