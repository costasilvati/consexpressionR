# Title     : Run baySeq
# Objective : Execute baySeq with default paraemters by a count table
# Created by: julianacostasilva
# Created on: 30/03/22

library("baySeq")
library("edgeR")
library("limma")

read_count_file <- function (table_count, separador){
    table <- read.csv( table_count,  sep = separador, row.names = 1, header = TRUE, stringsAsFactors = FALSE)
    # m <- as.matrix(table)
    return(table)
}

run_baySeq <- function  (m, replics, replicates, out){
    if(require("parallel")) cl <- makeCluster(4) else cl <- NULL
    groups <- list(NDE = replicates, DE = replicates)
    CD <- new("countData", data = m, replicates = replicates, groups = groups)
    libsizes(CD) <- getLibsizes(CD)
    CD <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl, equalDispersions = TRUE)
    CD <- getLikelihoods(CD, prs=c(0.5, 0.5), pET="BIC", cl=cl)
    write.table(topCounts(CD, group = "DE", number = 65000, normaliseData = TRUE), out, sep="\t", quote = FALSE)
    return (CD)
}

run_edger <- function (m, replics, replicates, out){
    group <- c(replicates)
    y.dge <- DGEList(counts = m, group = group)
    if (replics < 1){
        print('Replicates not found by edgeR. EdgeR should be executed manual form.')
    }else if(replics == 1){
        bcv <- 0.2
        y.et <- exactTest(y.dge, dispersion = bcv^2)
        y.tp <- topTags(y.et, n = 100000)
        # y.pvalues <- y.et$table$PValue
        write.table(y.tp$table, out, sep = "\t", quote = FALSE)
    }else{
        y.dge <- calcNormFactors(y.dge)
        y.dge <- estimateDisp(y.dge)
        y.dge <- estimateCommonDisp(y.dge)
        y.et <- exactTest(y.dge)
        y.tp <- topTags(y.et, n = 100000)
        y.pvalues <- y.et$table$PValue
        write.table(y.tp$table, out, sep = "\t", quote = FALSE)
    }
    return(y.tp$table)
}

run_limma <- function (m, replics, replicates, out, methodNorm = "TMM", methodAdjPvalue = "BH", numberTopTable = 1000000){
    if (replics <= 1){
        print('ERROR: limma-voom require more than one replics.')
    }else {
        nf <- calcNormFactors(m, method = methodNorm)
        condition = factor(c(replicates))
        voom.data <- voom(m, model.matrix(~factor(condition)), lib.ize = colSums(m) * nf)
        voom.data$genes = rownames(m)
        voom.fitlimma = lmFit(voom.data, design=model.matrix(~factor(condition)))
        voom.fitbayes = eBayes(voom.fitlimma)
        voom.pvalues = voom.fitbayes$p.value[, 2]
        voom.adjpvalues = p.adjust(voom.pvalues, method=methodAdjPvalue)
        design <- group
        data <- topTable(voom.fitbayes, coef=ncol(design), number=numberTopTable)
        write.table(data, file=out, sep = "\t", quote = FALSE)
        return(data)
    }
}


# Teste
replics <- 3
group_name <- c("0b","1b")
table_count = '/Users/julianacostasilva/Library/CloudStorage/OneDrive-Pessoal/ProjetoDoutorado/bioconvergencia/reads_RNApa/kallisto_quant_align_apa_1B_0B/gonoda/gonoda_tpm.csv'
separador = "\t"
m <- as.matrix(read_count_file(table_count, separador))
replicates <- rep(group_name, each = replics)
#out = '/Volumes/SD128/bioconvergencia/reads_RNApa/RNApa_apa_1B_0B_gonoda_consexpression_baySeq.csv'
# bayseq_de <- run_baySeq(m, replics, replicates, out)
# out = '/Volumes/SD128/bioconvergencia/reads_RNApa/RNApa_apa_1B_0B_gonoda_consexpression_edger.csv'
#edger_de <- run_edger(m,replics, replicates, out)
#out = '/Volumes/SD128/bioconvergencia/reads_RNApa/RNApa_apa_1B_0B_gonoda_consexpression_edger.csv'
limma_de <- run_limma(m,replics, replicates, out)
out = '/Volumes/SD128/bioconvergencia/reads_RNApa/RNApa_apa_1B_0B_gonoda_consexpression_limma.csv'
# install.packages("styler") # formata
# https://cran.r-project.org/doc/manuals/R-exts.html#Creating-R-packages
# http://bioconductor.org/developers/package-submission/
# https://www.bioconductor.org/developers/how-to/coding-style/


