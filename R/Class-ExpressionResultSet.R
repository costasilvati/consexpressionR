#' An S4 class to store differential expression results
#'
#' @slot results A named list of data.frames, each corresponding to a method (e.g., edgeR, DESeq2, etc.)
#' @slot methods A character vector of methods used
#' @slot parameters A list of parameters used in the analysis
#' @slot consensus A list containing consensus DEG results (e.g., from consensusList)
#'
#' @export
setClass("ExpressionResultSet",
         slots = list(
           results = "list",
           methods = "character",
           parameters = "list",
           consensus = "list"
         ))

#' Constructor for ExpressionResultSet
#'
#' @param results Named list of data.frames
#' @param methods Character vector with method names
#' @param parameters List of analysis parameters
#'
#' @return An ExpressionResultSet object
#' @export
createExpressionResultSet <- function(results, methods, parameters = list()) {
  # Verifica se todos os elementos são data.frames ou NULL
  if (!all(sapply(results, function(x) is.null(x) || inherits(x, "data.frame")))) {
    stop("All elements in 'results' must be data.frames or NULL.")
  }

  if (!is.character(methods)) {
    stop("'methods' must be a character vector.")
  }

  if (!is.list(parameters)) {
    stop("'parameters' must be a list.")
  }

  new("ExpressionResultSet",
      results = results,
      methods = methods,
      parameters = parameters)
}


#' Show method for ExpressionResultSet
#'
#' @param object An ExpressionResultSet object
#' @export
setMethod("show", "ExpressionResultSet", function(object) {
  cat("ExpressionResultSet object\n")
  cat("Methods executed:", paste(object@methods, collapse = ", "), "\n")
  cat("Results available:", length(object@results), "data.frames\n")

  if (length(object@consensus) == 0) {
    cat("Consensus: not computed or empty.\n")
  } else {
    cat("Consensus: available (", length(object@consensus), " gene(s) listed)\n")
  }
})

#' Summary method for ExpressionResultSet
#'
#' @param object An ExpressionResultSet object
#' @export
setMethod("summary", "ExpressionResultSet", function(object) {
  cat("Summary of ExpressionResultSet\n")
  cat("---------------------------------\n")
  cat("Methods used:", paste(object@methods, collapse = ", "), "\n")

  nTables <- length(object@results)
  cat("Number of result tables:", nTables, "\n")

  if (nTables > 0) {
    sizes <- sapply(object@results, function(x) if (is.null(x)) 0 else nrow(x))
    cat("Average number of DEGs per method:", round(mean(sizes), 2), "\n")
  }

  if (length(object@consensus) == 0) {
    cat("Consensus list is empty or not defined.\n")
  } else {
    cat("Number of genes in consensus:", length(object@consensus), "\n")
  }
})
