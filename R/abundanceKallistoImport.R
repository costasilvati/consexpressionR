#' Import kallisto abundance to DESeq2 data set
#'
#' @param dir path to kallisto results
#' @param pathReportRuns report from experiment contains columns pop : population,center: local from data,run: name for each run,condition: condition by run (ex.: treatment, treated))
#' @param pathReportFile parent directory that contains kallisto results organized by directory named like run
#' @param replics number of biological or technical replicates
#' @param dirRuns directory that contains kallisto results organized by directory named like run
#' @param fileKallistoAbundance abundabce.tsv
#' @param conditionList list by treatments per group
#'
#' @return DESeq2 DataSet
#' @export
#'
#' @examples
abundanceKallistoImport <- function(pathReportRuns,
                                    pathReportFile,
                                    dirRuns,
                                    fileKallistoAbundance,
                                    conditionList){
  count <- utils::read.table(file.path(pathReportRuns),
                             header=TRUE) # gonoda_report.txt
  count$condition <- conditionList
  rownames(count) <- count$run
  count[,c("pop","center","run","condition")]
  files <- file.path(pathReportFile,
                     dirRuns,
                     count$run,
                     fileKallistoAbundance)
  files
  names(files) <- count$run
  txi <- tximport::tximport(files,
                            type = "kallisto",
                            txOut = TRUE)
  sampleTable <- data.frame(condition = factor(count$run))
  rownames(sampleTable) <- colnames(txi$counts)
  ddsTxi <- DESeq2::DESeqDataSetFromTximport(txi,
                                             design = ~ 1,
                                             colData = sampleTable)
  return(ddsTxi)
}
