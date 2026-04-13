#' Make expression analysis of multiple tools and return results by tool.
#'
#' @importFrom utils read.csv write.csv
#' @importFrom methods new
#' @param numberReplics number of replicates (technical or biological) per group
#' @param groupName character vector, group/treatment names (e.g., c("BM","JJ"))
#' @param tableCountPath path to csv file that contains count/abundance data (local)
#' @param sepCharacter separator used to split csv data (comma or tab)
#' @param rDataFrameCount data.frame/matrix with table count; rownames are genes and colnames are samples
#' @param experimentName experiment name
#' @param outDirPath output directory (created if needed)
#' @param printResults if TRUE, writes per-tool outputs to outDirPath
#' @param methodNormLimma normalization method used by limma/edgeR
#' @param methodAdjPvalueLimma p-value adjustment method for limma
#' @param numberTopTableLimma max number of genes for limma topTable
#' @param methodNormEdgeR normalization method for edgeR
#' @param kNoiseq NOISeq parameter
#' @param factorNoiseq NOISeq parameter
#' @param lcNoiseq NOISeq parameter
#' @param replicatesNoiseq NOISeq parameter
#' @param normNoiseq NOISeq parameter
#' @param condExpNoiseq NOISeq parameter
#' @param respTypeSamseq SAMSeq parameter
#' @param npermSamseq SAMSeq parameter
#' @param fdrEbseq EBSeq parameter
#' @param maxRoundEbseq EBSeq parameter
#' @param methodDeResultsEbseq EBSeq parameter
#' @param filterIdKnowseq KnowSeq parameter
#' @param notSapiensKnowseq KnowSeq parameter
#' @param fitTypeDeseq2 DESeq2 fitType
#' @param controlDeseq2 control group label for DESeq2
#' @param contrastDeseq2 contrast group label for DESeq2
#' @param deNovoAanalysis if TRUE, skips KnowSeq (requires annotation)
#' @param progressShiny optional shiny progress callback
#'
#' @return An \code{ExpressionResultSet} object containing results per method.
#' @export
#'
#' @examples
#' data(gse95077)
#' treats <- c("BM", "JJ")
#' res <- runExpression(
#'   numberReplics = 3,
#'   groupName = treats,
#'   rDataFrameCount = gse95077,
#'   controlDeseq2 = "BM",
#'   contrastDeseq2 = "JJ",
#'   printResults = FALSE,
#'   outDirPath = tempdir()
#' )
#' summary(res)
runExpression <- function(numberReplics,
                          groupName,
                          tableCountPath = "data/gse95077.csv",
                          sepCharacter = ",",
                          rDataFrameCount = NULL,
                          experimentName = "genericExperiment",
                          outDirPath = tempdir(),
                          printResults = FALSE,
                          fitTypeDeseq2 = "local",
                          controlDeseq2 = "",
                          contrastDeseq2 = "",
                          methodNormLimma = "TMM",
                          methodAdjPvalueLimma = "BH",
                          numberTopTableLimma = 1000000,
                          filterIdKnowseq = "ensembl_gene_id",
                          notSapiensKnowseq = FALSE,
                          methodNormEdgeR = "TMM",
                          normNoiseq = "rpkm",
                          kNoiseq = 0.5,
                          factorNoiseq = "Tissue",
                          lcNoiseq = 0,
                          replicatesNoiseq = "technical",
                          condExpNoiseq = c(""),
                          respTypeSamseq = "Two class unpaired",
                          npermSamseq = 100,
                          fdrEbseq = 0.05,
                          maxRoundEbseq = 50,
                          methodDeResultsEbseq = "robust",
                          deNovoAanalysis = FALSE,
                          progressShiny = NULL) {

    .progress <- function(detail) {
        if (!is.null(progressShiny)) {
            progressShiny(detail = detail)
        }
    }
    .is_count_like <- function(x) {
        x <- as.matrix(x)
        if (anyNA(x) || any(x < 0)) {
            return(FALSE)
        }
        all(abs(x - round(x)) < .Machine$double.eps^0.5)
    }
    if (!is.null(rDataFrameCount)) {
        countMatrix <- as.matrix(rDataFrameCount)
    } else {
        countMatrix <- as.matrix(readCountFile(tableCountPath, sepCharacter))
    }
    .progress("Count matrix")

    if (is.null(colnames(countMatrix))) {
        stop("countMatrix must have column names (sample names).")
    }
    if (!is.numeric(numberReplics) || length(numberReplics) != 1L || numberReplics < 1) {
        stop("'numberReplics' must be a single number >= 1.")
    }
    designExperiment <- rep(groupName, each = numberReplics)
    if (length(designExperiment) != ncol(countMatrix)) {
        stop(sprintf(
            "Expected %d samples from groupName/numberReplics but countMatrix has %d columns.",
            length(designExperiment), ncol(countMatrix)
        ))
    }
    # ---- run tools ----
    resultTool <- list()
    .progress("Executing edgeR...")
    resultTool$edger <- tryCatch(
        runEdger(countMatrix, numberReplics, designExperiment, methNorm = methodNormEdgeR),
        error = function(e) {
            warning(sprintf("edgeR failed and will be skipped: %s", conditionMessage(e)))
            NULL
        }
    )
    .progress("Executing KnowSeq...")
    resultTool$knowseq <- NULL
    if (deNovoAanalysis) {
        warning("KnowSeq skipped for de novo analysis (annotation not available).")
    } else {
        resultTool$knowseq <- tryCatch(
            runKnowSeq(
                count = countMatrix,
                groupName = groupName,
                numberReplic = numberReplics,
                filterId = filterIdKnowseq,
                notSapiens = notSapiensKnowseq
            ),
            error = function(e) {
                warning(sprintf("KnowSeq failed and will be skipped: %s", conditionMessage(e)))
                NULL
            }
        )
    }
    .progress("Executing limma...")
    resultTool$limma <- tryCatch(
        runLimma(
            countMatrix, numberReplics, designExperiment,
            methodNormLimma, methodAdjPvalueLimma, numberTopTableLimma
        ),
        error = function(e) {
            warning(sprintf("limma failed and will be skipped: %s", conditionMessage(e)))
            NULL
        }
    )
    .progress("Executing NOISeq...")
    resultTool$noiseq <- tryCatch(
        runNoiSeq(
            countMatrix = countMatrix,
            groups = groupName,
            designExperiment = designExperiment,
            normParm = normNoiseq,
            kParam = kNoiseq,
            factorParam = factorNoiseq,
            lcParam = lcNoiseq,
            replicatesParam = replicatesNoiseq,
            condExp = condExpNoiseq
        ),
        error = function(e) {
            warning(sprintf("NOISeq failed and will be skipped: %s", conditionMessage(e)))
            NULL
        }
    )
    .progress("Executing EBSeq...")
    resultTool$ebseq <- tryCatch(
        runEbseq(
            as.matrix(countMatrix),
            designExperiment,
            fdr = fdrEbseq,
            maxRound = maxRoundEbseq,
            methodDeResults = methodDeResultsEbseq,
            groups = groupName
        ),
        error = function(e) {
            warning(sprintf("EBSeq failed and will be skipped: %s", conditionMessage(e)))
            NULL
        }
    )
    if (!.is_count_like(countMatrix)) {
        warning("DESeq2 and SAMSeq skipped: input does not look like raw count data (non-integer values detected).")
        resultTool$deseq2 <- NULL
        resultTool$samseq <- NULL
    } else {
        .progress("Executing DESeq2...")
        resultTool$deseq2 <- tryCatch(
            runDeseq2(
                countMatrix,
                groupName,
                numberReplics,
                controlGroup = controlDeseq2,
                contrastGroup = contrastDeseq2,
                fitTypeParam = fitTypeDeseq2
            ),
            error = function(e) {
                warning(sprintf("DESeq2 failed and will be skipped: %s", conditionMessage(e)))
                NULL
            }
        )
        .progress("Executing SAMSeq...")
        resultTool$samseq <- tryCatch(
            runSamSeq(
                countMatrix,
                designExperiment,
                respType = respTypeSamseq,
                numberPermutations = npermSamseq
            ),
            error = function(e) {
                warning(sprintf("SAMSeq failed and will be skipped: %s", conditionMessage(e)))
                NULL
            }
        )
    }
    if (isTRUE(printResults)) {
        .progress("Writing results...")
        if (!dir.exists(outDirPath)) {
            dir.create(outDirPath, recursive = TRUE, showWarnings = FALSE)
        }
        tryCatch(
            {
                for (nm in names(resultTool)) {
                    data <- resultTool[[nm]]
                    if (is.null(data)) next
                    utils::write.csv(
                        data,
                        file = file.path(outDirPath, paste0(experimentName, "_", nm, ".csv"))
                    )
                    save(
                        data,
                        file = file.path(outDirPath, paste0(experimentName, "_", nm, ".RData"))
                    )
                }
                toolsFramed <- frameAllGenes(resultTool, countMatrix)
                save(
                    toolsFramed,
                    file = file.path(outDirPath, paste0(experimentName, "_toolsResult.RData"))
                )
            },
            error = function(e) {
                warning(sprintf("Writing results failed: %s", conditionMessage(e)))
            }
        )
    }
    .progress("Complete!")
    message("Differential expression analysis completed.")
    parameters <- list(
        numberReplics = numberReplics,
        groupName = groupName,
        experimentName = experimentName,
        tableCountPath = tableCountPath,
        printResults = printResults
    )
    createExpressionResultSet(
        results = resultTool,
        methodNames = names(resultTool),
        parameters = parameters,
        consensus = list()
    )
}
