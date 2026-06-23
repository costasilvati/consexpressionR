#'  Execute SAMSeq gene Expression analysis to count data
#'
#' @param countMatrix either a matrix of raw (read) counts.
#' @param respType Problem type: "Quantitative" for a continuous parameter; "Two class unpaired" for two classes with unpaired observations; "Survival" for censored survival outcome; "Multiclass": more than 2 groups; "Two class paired" for two classes with paired observations.
#' @param designExperiment replicate and treatment by samples
#' @param numberPermutations Number of permutations used to estimate false discovery rates
#'
#' @return SAMSeq report in data Frame
#' @export
#'
#' @examples
#' set.seed(42)
#' counts <- matrix(
#'   as.integer(c(
#'     rnbinom(200, mu = 10,  size = 1), 
#'     rnbinom(200, mu = 100, size = 1) 
#'   )),
#'   nrow = 100,
#'   dimnames = list(paste0("gene", seq_len(100)),
#'                     paste0("sample", seq_len(4)))
#' )
#' groups_info <- c("control", "control", "treated", "treated")
#' treats <- c("control", "treated")
#' toolResult <- NULL
#' if (requireNamespace("samr", quietly = TRUE)) {
#'   toolResult$samseq <- runSamSeq(countMatrix = counts,
#'                                designExperiment = rep(treats, each = 2))
#' }
runSamSeq <- function (countMatrix, designExperiment, respType="Two class unpaired", numberPermutations=100){
    .check_package("samr", repo = "CRAN")
    samResult <- samr::SAMseq(countMatrix, as.factor(designExperiment), resp.type=respType,
                              geneid=row.names(countMatrix), genenames=row.names(countMatrix), nperms=numberPermutations)
    samResultTable <- rbind(samResult$siggenes.table$genes.up, samResult$siggenes.table$genes.lo)
    samScore <- rep(0,nrow(countMatrix))
    samScore[match(samResultTable[,1], rownames(countMatrix))] <- as.numeric(samResultTable[,3])
    samFdr <- rep(1, nrow(countMatrix))
    samFdr[match(samResultTable[,1], rownames(countMatrix))] <- as.numeric(samResultTable[,5])/100
    rownames(samResultTable) <- samResultTable[,1]
    result<- as.data.frame(samResultTable)
    return(result)
}
