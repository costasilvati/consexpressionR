# Title     : run DESeq
# Objective : execute DESeq by informations user
# Created by: julianacostasilva
# Created on: 10/03/22

# Fonte: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#tximport

library("tximport")
library("readr")

BiocManager::install("tximportData")
library("tximportData")
dir <- system.file("extdata", package="tximportData")
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
samples$condition <- factor(rep(c("A","B"),each=3))
rownames(samples) <- samples$run
samples[,c("pop","center","run","condition")]

