# Title     : Run baySeq
# Objective : Execute baySeq with default paraemters by a count table
# Created by: julianacostasilva
# Created on: 30/03/22

# install.packages("styler") # formata
# https://cran.r-project.org/doc/manuals/R-exts.html#Creating-R-packages
# http://bioconductor.org/developers/package-submission/
# https://www.bioconductor.org/developers/how-to/coding-style/
# Livro R Packages - https://r-pkgs.org/




consexpression2 <- function (numberReplics=3, 
                             groupName= c("0b", "1b"),
                             tableCountPath="table_count.csv",
                             split="\t",
                             experimentName="genericExperiment",
                             outDirPath="."){
    # validações?
    countMatrix <- readCountFile(tableCountPath, split)
    designExperiment <- rep(groupName, each = numberReplics)
    # baySeq ok
    #bayseqResult<-runBaySeq(countMatrix, replicates, createNameFileOutput(outDirPath,experimentName,execName='baySeq'))
    # edgeR ok
    #edgerResult<-runEdger(countMatrix, numberReplics, designExperiment, createNameFileOutput(outDirPath,experimentName,execName='edgeR'))
    # limma ok
    #limmaResult<-runLimma(countMatrix, numberReplics, designExperiment, createNameFileOutput(outDirPath,experimentName,execName='limma'))
    # NOISeq ok
    #noiseqResult<-runNOISeq(countMatrix, designExperiment, createNameFileOutput(outDirPath,experimentName,execName='NOISeq'))
    #samseqResult<-runSamSeq(countMatrix, numberReplics, designExperiment, createNameFileOutput(outDirPath,experimentName,execName='SAMSeq'))
    # EBSeq ??s
    # ebseqResult <- runEbseq(countMatrix, numberReplics, designExperiment, createNameFileOutput(outDirPath,experimentName,execName = 'EBSeq'))
    # DESeq ??
    # deseqResult <- runDeseq(countMatrix, numberReplics, designExperiment, createNameFileOutput(outDirPath,experimentName,execName = 'DESeq'))
    # DESeq2 ??
    # deseq2OutPath = '/Volumes/SD128/bioconvergencia/reads_RNApa/RNApa_apa_1B_0B_gonoda_consexpression_deseq2.csv'
    # deseq2Result <- runDeseq(countMatrix, numberReplics, designExperiment, deseq2OutPath)
}
