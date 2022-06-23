# pathReportRuns<- "gonoda_report.txt"
# pathReportFile <-"/Users/julianacostasilva/Library/CloudStorage/OneDrive-Pessoal/ProjetoDoutorado/bioconvergencia/reads_RNApa"
# conditionList <- c("0b","1b")
# replics <- 3

abundanceKallistoImport <- function(pathReporRuns,pathReportFile,conditionList, dirAbundanceFiles, replics){
  count <- read.table(file.path(pathReportRuns), header=TRUE)
  count$condition <- factor(rep(conditionList,each=replics))
  rownames(count) <- count$run
  count[,c("pop","center","run","condition")]
  files <- file.path(pathReportFile, "kallisto", count$run, "abundance.tsv")
  names(files) <- paste0("sample", 1:6)
  txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
  head(txi.kallisto.tsv$counts)

}
