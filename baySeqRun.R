# Title     : Run baySeq
# Objective : Execute baySeq with default paraemters by a count table
# Created by: julianacostasilva
# Created on: 30/03/22

library("baySeq")

run_baySeq <- function  (table_count, separador, replics, group_name, output){
   if(require("parallel")) cl <- makeCluster(4) else cl <- NULL
   replicates <- rep(group_name, each = replics)
   table <- read.csv( table_count,  sep = separador, row.names = 1, header = TRUE, stringsAsFactors = FALSE)
   m <- as.matrix(table)
   groups <- list(NDE = c(1, (length(group_name) -1)))
   CD <- new("countData", data = m, replicates = replicates, groups = groups)
   libsizes(CD) <- getLibsizes(CD)
   CD <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl, equalDispersions = TRUE)
   CD <- getLikelihoods(CD, prs=c(0.5, 0.5), pET="BIC", cl=cl)
   write.table(topCounts(CD, group = "DE", number = 65000, normaliseData = TRUE), output, sep="\t", quote = FALSE)
}

# Teste
replics <- 3
group_name <- c("0b","1b")
column_name = c("01gon","05gon","09gon","04gon","04mus","05mus","06gon","06mus","08gon","08mus","09mus")
table_count = '/Volumes/SD128/bioconvergencia/reads_RNApa/kallisto_quant_RNApa_asca_0B_pb_tpm_tabular.csv'
separador = ","
out = '/Volumes/SD128/bioconvergencia/reads_RNApa/RNApa_apa_1B_0B-consexpression_baySeq.csv'
run_baySeq(table_count,separador, replics, group_name, out)


