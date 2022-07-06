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
    samResult <- samr::SAMseq(countMatrix,as.factor(designExperiment),resp.type='Multiclass',
                          geneid=row.names(countMatrix),genenames=row.names(countMatrix),nperms=100)
    samResultTable <- rbind(samResult$siggenes.table$genes.up, samResult$siggenes.table$genes.lo)
    samScore <- rep(0, nrow(countMatrix))
    samScore[match(samResultTable[,1], rownames(countMatrix))]=as.numeric(samResultTable[,3])
    samFdr = rep(1, nrow(countMatrix))
    samFdr[match(samResultTable[,1], rownames(countMatrix))] = as.numeric(samResultTable[,5])/100
    utils::write.table(samResultTable, file=samseqOutPtah, sep='\t', quote=FALSE)

    return(samResultTable)
}
