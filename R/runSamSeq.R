#'  Execute SAMSeq gene Expression analysis to count data
#'
#' @param countMatrix either a matrix of raw (read) counts.
#' @param replic number of replicate (technical or biologcal) integer
#' @param desingExperiment replicate and treatment by samples
#' @param samseqOutPath path to write SAMSeq report in csv file format
#'
#' @return SAMSeq report in data Frame
#' @export
#'
#' @examples
runSamSeq <- function (countMatrix, designExperiment, replic, samseqOutPtah){
    # SAMseq.test <- SAMseq(countMatrix, designExperiment, resp.type="Multiclass", nperms = 100)
    SAMseq.test <- SAMseq(countMatrix,as.factor(designExperiment),resp.type='Multiclass',
                          geneid=row.names(countMatrix),genenames=row.names(countMatrix),nperms=100)
    SAMseq.result.table <- rbind(SAMseq.test$siggenes.table$genes.up, SAMseq.test$siggenes.table$genes.lo)
    SAMseq.score <- rep(0, nrow(countMatrix))
    SAMseq.score[match(SAMseq.result.table[,1], rownames(countMatrix))]=as.numeric(SAMseq.result.table[,3])
    SAMseq.FDR = rep(1, nrow(countMatrix))
    SAMseq.FDR[match(SAMseq.result.table[,1], rownames(countMatrix))] = as.numeric(SAMseq.result.table[,5])/100
    write.table(SAMseq.result.table, file=samseqOutPtah, sep='\t', quote=FALSE)
}
