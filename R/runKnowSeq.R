#' Run differential expression analysis with KnowSeq
#'
#' Executes a two-group differential expression analysis using functions from
#' the \pkg{KnowSeq} package.
#'
#' @param count A matrix or data.frame of raw counts or abundance values, with
#'   genes in rows and samples in columns.
#' @param groupName A character vector with exactly two group labels.
#' @param numberReplic A single positive integer indicating the number of
#'   replicates per group.
#' @param filterId A single character string specifying the annotation filter
#'   used by \code{KnowSeq::getGenesAnnotation()}.
#' @param notSapiens A single logical value indicating whether a non-human
#'   annotation dataset should be used.
#'
#' @details
#' This function wraps the main \pkg{KnowSeq} workflow steps:
#' \enumerate{
#'   \item retrieve gene annotation with
#'     \code{KnowSeq::getGenesAnnotation()};
#'   \item compute expression values with
#'     \code{KnowSeq::calculateGeneExpressionValues()};
#'   \item extract differentially expressed genes with
#'     \code{KnowSeq::DEGsExtraction()}.
#' }
#'
#' The function expects exactly two groups. Row names in \code{count} must
#' contain gene identifiers compatible with the selected \code{filterId}.
#'
#' Because this workflow depends on external package behavior and may require
#' internet access, errors raised by internal \pkg{KnowSeq}, \pkg{cqn}, or
#' \pkg{mclust} calls are rethrown with additional context.
#'
#' @return A data.frame containing the differential expression results returned
#'   by \code{KnowSeq::DEGsExtraction()}.
#'
#' @seealso
#' \code{\link{runExpression}},
#' \code{\link{expressionDefinition}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(gse95077)
#' treats <- c("BM", "JJ")
#' res <- runKnowSeq(
#'   count = gse95077,
#'   groupName = treats,
#'   numberReplic = 3
#' )
#' head(res)
#' }
runKnowSeq <- function(
        count,
        groupName,
        numberReplic,
        filterId = "ensembl_gene_id",
        notSapiens = FALSE
) {
    if (!is.matrix(count) && !is.data.frame(count)) {
        stop(
            sprintf(
                "'count' must be a matrix or data.frame; got an object of class %s.",
                paste(class(count), collapse = ", ")
            ),
            call. = FALSE
        )
    }

    if (!is.character(groupName)) {
        stop(
            sprintf(
                "'groupName' must be a character vector; got an object of class %s.",
                paste(class(groupName), collapse = ", ")
            ),
            call. = FALSE
        )
    }

    if (length(groupName) != 2L) {
        stop(
            sprintf(
                "'groupName' must contain exactly 2 group labels; got %d value(s).",
                length(groupName)
            ),
            call. = FALSE
        )
    }

    if (!is.numeric(numberReplic) || length(numberReplic) != 1L || is.na(numberReplic)) {
        stop(
            sprintf(
                "'numberReplic' must be a single non-missing numeric value; got length %d.",
                length(numberReplic)
            ),
            call. = FALSE
        )
    }

    if (numberReplic < 1) {
        stop(
            sprintf(
                "'numberReplic' must be >= 1; got %s.",
                format(numberReplic)
            ),
            call. = FALSE
        )
    }

    if (!isTRUE(all.equal(numberReplic, as.integer(numberReplic)))) {
        stop(
            sprintf(
                "'numberReplic' must be an integer value; got %s.",
                format(numberReplic)
            ),
            call. = FALSE
        )
    }

    if (!is.character(filterId) || length(filterId) != 1L || is.na(filterId)) {
        stop(
            sprintf(
                "'filterId' must be a single non-missing character string; got length %d.",
                length(filterId)
            ),
            call. = FALSE
        )
    }

    if (!is.logical(notSapiens) || length(notSapiens) != 1L || is.na(notSapiens)) {
        stop(
            sprintf(
                "'notSapiens' must be a single non-missing logical value; got length %d.",
                length(notSapiens)
            ),
            call. = FALSE
        )
    }

    required_pkgs <- c("KnowSeq", "cqn", "mclust")
    pkg_available <- vapply(
        required_pkgs,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
    missing_pkgs <- required_pkgs[!pkg_available]

    if (length(missing_pkgs) > 0L) {
        stop(
            sprintf(
                "Missing required package(s): %s.",
                paste(missing_pkgs, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    count_matrix <- as.matrix(count)
    gene_names <- rownames(count_matrix)

    if (is.null(gene_names)) {
        stop(
            sprintf(
                "'count' must have row names containing gene identifiers; no row names were found."
            ),
            call. = FALSE
        )
    }

    if (anyNA(gene_names)) {
        stop(
            sprintf(
                "'count' row names must not contain missing values; detected %d missing value(s).",
                sum(is.na(gene_names))
            ),
            call. = FALSE
        )
    }

    if (any(gene_names == "")) {
        stop(
            sprintf(
                "'count' row names must not contain empty strings; detected %d empty value(s).",
                sum(gene_names == "")
            ),
            call. = FALSE
        )
    }

    numberReplic <- as.integer(numberReplic)
    design_experiment <- rep(groupName, each = numberReplic)
    expected_samples <- length(design_experiment)
    observed_samples <- ncol(count_matrix)

    if (observed_samples != expected_samples) {
        stop(
            sprintf(
                "'count' has %d column(s), but 'groupName' and 'numberReplic' imply %d sample(s).",
                observed_samples,
                expected_samples
            ),
            call. = FALSE
        )
    }

    result <- tryCatch(
        {
            annotation <- KnowSeq::getGenesAnnotation(
                values = as.character(gene_names),
                filter = filterId,
                notHSapiens = notSapiens
            )

            expression_matrix <- KnowSeq::calculateGeneExpressionValues(
                countsMatrix = count_matrix,
                annotation = annotation,
                genesNames = FALSE
            )

            knowseq_result <- KnowSeq::DEGsExtraction(
                expression_matrix,
                labels = design_experiment,
                multiDegsMethod = "cov",
                lfc = -20
            )

            if (!is.list(knowseq_result)) {
                stop(
                    sprintf(
                        "KnowSeq returned an object of class %s; expected a list-like result.",
                        paste(class(knowseq_result), collapse = ", ")
                    ),
                    call. = FALSE
                )
            }

            if (!"DEG_Results" %in% names(knowseq_result)) {
                stop(
                    sprintf(
                        "KnowSeq returned an unexpected result: component 'DEG_Results' was not found."
                    ),
                    call. = FALSE
                )
            }

            if (!"DEGs_Table" %in% names(knowseq_result$DEG_Results)) {
                stop(
                    sprintf(
                        "KnowSeq returned an unexpected result: component 'DEG_Results$DEGs_Table' was not found."
                    ),
                    call. = FALSE
                )
            }

            knowseq_result$DEG_Results$DEGs_Table
        },
        error = function(e) {
            msg <- conditionMessage(e)

            if (grepl("mclustBIC|Mclust", msg)) {
                stop(
                    sprintf(
                        "KnowSeq execution failed due to an internal dependency issue involving the 'mclust' package: %s",
                        msg
                    ),
                    call. = FALSE
                )
            }

            stop(
                sprintf(
                    "runKnowSeq() failed during KnowSeq execution: %s",
                    msg
                ),
                call. = FALSE
            )
        }
    )

    result
}
