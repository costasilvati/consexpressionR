#' Run differential expression analysis using baySeq
#'
#' Performs differential expression analysis on RNA-seq count data using
#' the baySeq package, which applies empirical Bayesian methods to estimate
#' posterior likelihoods of differential expression.
#'
#' @param data A numeric matrix of raw read counts, with genes as rows and
#'   samples as columns. Row names should contain gene identifiers.
#' @param groups A named list defining the group structure for baySeq.
#'   For a two-group comparison, use:
#'   \code{list(NDE = c(1,1,1,1), DE = c(1,1,2,2))}, where values indicate
#'   group membership for each sample. If \code{NULL}, groups are inferred
#'   from \code{sampleInfo}.
#' @param sampleInfo A character or factor vector indicating group membership
#'   for each sample (e.g., \code{c("control","control","treated","treated")}).
#'   Used to build \code{groups} automatically when \code{groups = NULL}.
#' @param pCutoff Numeric. Posterior probability cutoff to call a gene as
#'   differentially expressed. Default is \code{0.95}.
#' @param fdrCutoff Numeric. FDR cutoff applied to adjusted p-values for
#'   secondary filtering. Default is \code{0.05}.
#' @param samplesize Integer. Number of samples used to estimate priors in
#'   baySeq. Default is \code{10000}. Reduce for faster execution during
#'   testing.
#' @param seed Integer. Random seed for reproducibility. Default is \code{42}.
#'
#' @return A data frame with one row per gene, containing:
#'   \describe{
#'     \item{gene}{Gene identifier (from row names of \code{data})}
#'     \item{posterior_DE}{Posterior probability of differential expression}
#'     \item{FDR}{False discovery rate estimate}
#'     \item{logFC}{Log2 fold-change between groups
#'       (mean group 2 / mean group 1)}
#'     \item{DE_baySeq}{Logical. \code{TRUE} if the gene is called as
#'       differentially expressed based on \code{pCutoff} and
#'       \code{fdrCutoff}}
#'   }
#'   Returns \code{NULL} invisibly if baySeq is not installed or if the
#'   analysis fails.
#'
#' @details
#' The function wraps the baySeq workflow for two-group RNA-seq comparisons.
#' It constructs a \code{countData} object, estimates prior distributions via
#' \code{getPriors.NB}, and computes likelihoods via \code{getLikelihoods}.
#' Genes are called as differentially expressed when their posterior
#' probability of DE exceeds \code{pCutoff} AND their FDR is below
#' \code{fdrCutoff}.
#'
#' The log2 fold-change is calculated from the group means of the raw count
#' matrix and is provided for reference only; baySeq itself does not model
#' fold-change directly.
#'
#' @references
#' Hardcastle TJ and Kelly KA (2010). baySeq: Empirical Bayesian methods
#' for identifying differential expression in sequence count data.
#' BMC Bioinformatics, 11, 422. \doi{10.1186/1471-2105-11-422}
#'
#' @seealso \code{\link[baySeq]{countData}}, \code{\link[baySeq]{getPriors.NB}},
#'   \code{\link[baySeq]{getLikelihoods}}
#'
#' @importFrom stats aggregate
#'
#' @examples
#' if (requireNamespace("baySeq", quietly = TRUE)) {
#'   set.seed(42)
#'   counts <- matrix(
#'     as.integer(c(
#'       rnbinom(200, mu = 10,  size = 1), 
#'       rnbinom(200, mu = 100, size = 1) 
#'     )),
#'     nrow = 100,
#'     dimnames = list(paste0("gene", seq_len(100)),
#'                     paste0("sample", seq_len(4)))
#'   )
#'   groups_info <- c("control", "control", "treated", "treated")
#'   result <- runBayseq(counts, sampleInfo = groups_info)
#'   utils::head(result)
#' }
#'
#' @export
runBayseq <- function(data,
                      groups = NULL,
                      sampleInfo  = NULL,
                      pCutoff     = 0.95,
                      fdrCutoff   = 0.05,
                      samplesize  = 10000,
                      seed        = 42) {
    .check_package("baySeq", repo = "Bioconductor")
    if (!is.matrix(data)) {
        data <- as.matrix(data)
    }
    if (!is.integer(data)) {
        message("[runBayseq] Converting matrix to integer values (this is necessary to execute baySeq).")
        mode(data) <- "integer"
    }
    if (is.null(rownames(data))) {
        rownames(data) <- paste0("gene", seq_len(nrow(data)))
    }
    if (is.null(groups)) {
        if (is.null(sampleInfo)) {
            stop("[runBayseq] Please, infrom 'groups' or 'sampleInfo' to define groups.")
        }
        group_factor <- as.factor(sampleInfo)
        group_levels <- levels(group_factor)
        if (length(group_levels) != 2) {
            stop("[runBayseq] consexpressionR has suport only two groups compair.")
        }
        group_vec <- as.factor(as.integer(group_factor))
        groups <- list(
            NDE = rep(1L, ncol(data)),
            DE  = group_vec
        )
    }
    replicates_factor <- factor(groups$DE)
    
    if (length(levels(replicates_factor)) != 2) {
        stop("[runBayseq] 'groups$DE' must define exactly two groups.")
    }
    cd <- methods::new(
        "countData",
        data   = data,
        replicates = factor(sampleInfo),
        groups = groups
    )
    message("baySeq libsizes")
    baySeq::libsizes(cd) <- baySeq::getLibsizes(cd)
    set.seed(seed)
    message("baySeq getPriors")
    cd <- tryCatch(
        baySeq::getPriors.NB(cd, samplesize = samplesize, estimation = "QL",
                             cl = NULL),
        error = function(e) {
            warning("[runBayseq] Fail in priors estimation: ", conditionMessage(e))
            return(NULL)
        }
    )
    
    if (is.null(cd)) return(invisible(NULL))
    cd <- tryCatch(
        baySeq::getLikelihoods(cd, pET = "BIC", cl = NULL),
        error = function(e) {
            warning("[runBayseq] Fail in calc likelihoods: ", conditionMessage(e))
            return(NULL)
        }
    )
    if (is.null(cd)) return(invisible(NULL))
    results_raw <- baySeq::topCounts(cd, group = "DE", number = nrow(data),
                                     normaliseData = TRUE)
    cat("Available columns in results_raw:\n")
    print(colnames(results_raw))
    print(utils::head(results_raw))
    grp1_idx <- which(groups$DE == 1L)
    grp2_idx <- which(groups$DE == 2L)
    mean_grp1 <- rowMeans(data[, grp1_idx, drop = FALSE] + 0.5)
    mean_grp2 <- rowMeans(data[, grp2_idx, drop = FALSE] + 0.5)
    logFC_vec  <- log2(mean_grp2 / mean_grp1)
    gene_names <- rownames(results_raw)
    result_df <- data.frame(
        gene         = gene_names,
        posterior_DE = results_raw[, "DE"],
        FDR          = results_raw[, "FDR.DE"],
        logFC        = logFC_vec[gene_names],
        DE_baySeq    = (results_raw[, "DE"] >= pCutoff) &
            (results_raw[, "FDR.DE"] <= fdrCutoff),
        row.names    = NULL,
        stringsAsFactors = FALSE
    )
    return(result_df)
}