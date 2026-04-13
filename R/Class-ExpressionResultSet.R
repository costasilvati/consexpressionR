#' An S4 class to store differential expression results
#'
#' @slot results A named list of data.frames, each corresponding to a method
#'   (e.g., edgeR, DESeq2, etc.)
#' @slot methodNames A character vector of methods used
#' @slot parameters A list of parameters used in the analysis
#' @slot consensus A list containing consensus DEG results
#'
#' @importClassesFrom methods classRepresentation
#' @importFrom methods new setClass setMethod setGeneric setValidity validObject show
#' @exportClass ExpressionResultSet
#' @examples
#' obj <- createExpressionResultSet(
#'   results = list(edgeR = data.frame(gene = "g1", logFC = 1)),
#'   methodNames = "edgeR",
#'   parameters = list(),
#'   consensus = list()
#' )
#' obj

setClass(
    "ExpressionResultSet",
    slots = list(
        results = "list",
        methodNames = "character",
        parameters = "list",
        consensus = "list"
    )
)

# ---- validity ----
#' @importFrom methods setValidity
setValidity("ExpressionResultSet", function(object) {
    # results: list, ideally named
    if (!is.list(object@results)) {
        return("'results' must be a list.")
    }

    # every element must be NULL or data.frame
    ok <- all(vapply(
        object@results,
        function(x) is.null(x) || inherits(x, "data.frame"),
        logical(1)
    ))
    if (!ok) {
        return("All elements in 'results' must be data.frames or NULL.")
    }

    # methodNames
    if (!is.character(object@methodNames)) {
        return("'methodNames' must be a character vector.")
    }
    if (length(object@methodNames) > 0L && anyNA(object@methodNames)) {
        return("'methodNames' cannot contain NA.")
    }

    # parameters
    if (!is.list(object@parameters)) {
        return("'parameters' must be a list.")
    }

    # consensus
    if (!is.list(object@consensus)) {
        return("'consensus' must be a list.")
    }

    TRUE
})

# ---- accessors ----
#' Get result tables
#'
#' @param object An ExpressionResultSet object.
#' @return A named list of data.frames (or NULL elements).
#' @export
#' @examples
#' obj <- createExpressionResultSet(
#'   results = list(edgeR = data.frame(gene = "g1", logFC = 1)),
#'   methodNames = "edgeR"
#' )
#' results(obj)

setGeneric("results", function(object) standardGeneric("results"))

#' @export
setMethod("results", "ExpressionResultSet", function(object) object@results)

#' Get method names used in the analysis
#'
#' @param object An ExpressionResultSet object.
#' @return A character vector.
#' @export
#' @examples
#' obj <- createExpressionResultSet(
#'   results = list(edgeR = data.frame(gene = "g1", logFC = 1)),
#'   methodNames = c("edgeR")
#' )
#' methodNames(obj)

setGeneric("methodNames", function(object) standardGeneric("methodNames"))

#' @export
setMethod("methodNames", "ExpressionResultSet", function(object) object@methodNames)

#' Get parameters used in the analysis
#'
#' @param object An ExpressionResultSet object.
#' @return A list.
#' @export
#' @examples
#' obj <- createExpressionResultSet(
#'   results = list(edgeR = data.frame(gene = "g1", logFC = 1)),
#'   methodNames = "edgeR",
#'   parameters = list(alpha = 0.01)
#' )
#' parameters(obj)

setGeneric("parameters", function(object) standardGeneric("parameters"))

#' @export
setMethod("parameters", "ExpressionResultSet", function(object) object@parameters)

#' Get consensus results
#'
#' @param object An ExpressionResultSet object.
#' @return A list.
#' @export
#' @examples
#' obj <- createExpressionResultSet(
#'   results = list(edgeR = data.frame(gene = c("g1","g2"), logFC = c(1,-1))),
#'   methodNames = "edgeR",
#'   consensus = list(g1 = TRUE)
#' )
#' consensus(obj)

setGeneric("consensus", function(object) standardGeneric("consensus"))

#' @export
setMethod("consensus", "ExpressionResultSet", function(object) object@consensus)

# ---- constructor ----
#' Constructor for ExpressionResultSet
#'
#' @param results Named list of data.frames (or NULL elements).
#' @param methodNames Character vector with method names.
#' @param parameters List of analysis parameters.
#' @param consensus List with consensus results (default: empty list).
#'
#' @return An ExpressionResultSet object.
#' @export
#' @examples
#' res_list <- list(
#'   edgeR = data.frame(gene = c("g1","g2"), logFC = c(1.2, -0.4)),
#'   DESeq2 = NULL
#' )
#'
#' obj <- createExpressionResultSet(
#'   results = res_list,
#'   methodNames = c("edgeR", "DESeq2"),
#'   parameters = list(alpha = 0.05),
#'   consensus = list(g1 = TRUE)
#' )
#'
#' obj
#' summary(obj)

createExpressionResultSet <- function(results,
                                      methodNames,
                                      parameters = list(),
                                      consensus = list()) {
    if (!is.list(results)) {
        stop("'results' must be a list.")
    }

    if (!all(vapply(results,
                    function(x) is.null(x) || inherits(x, "data.frame"),
                    logical(1)))) {
        stop("All elements in 'results' must be data.frames or NULL.")
    }

    if (!is.character(methodNames)) {
        stop("'methodNames' must be a character vector.")
    }

    if (!is.list(parameters)) {
        stop("'parameters' must be a list.")
    }

    if (!is.list(consensus)) {
        stop("'consensus' must be a list.")
    }

    obj <- new(
        "ExpressionResultSet",
        results = results,
        methodNames = methodNames,
        parameters = parameters,
        consensus = consensus
    )

    validObject(obj)
    obj
}

# ---- show ----
#' @importFrom methods setMethod show

setMethod("show", "ExpressionResultSet", function(object) {
    message("ExpressionResultSet object")

    m <- methodNames(object)
    message(sprintf("Methods executed: %s",
                    if (length(m) == 0L) "<none>" else paste(m, collapse = ", ")))

    res_list <- results(object)
    message(sprintf("Result tables: %d", length(res_list)))

    cons <- consensus(object)
    if (length(cons) == 0L) {
        message("Consensus: not computed or empty.")
    } else {
        message(sprintf("Consensus: available (%d gene(s))", length(cons)))
    }
})

# ---- summary ----

setMethod("summary", "ExpressionResultSet", function(object) {
    message("Summary of ExpressionResultSet")
    message("---------------------------------")

    m <- methodNames(object)
    message(sprintf("Methods used: %s",
                    if (length(m) == 0L) "<none>" else paste(m, collapse = ", ")))

    res_list <- results(object)
    nTables <- length(res_list)
    message(sprintf("Number of result tables: %d", nTables))

    if (nTables > 0L) {
        sizes <- vapply(
            res_list,
            function(x) if (is.null(x)) 0L else nrow(x),
            integer(1)
        )
        message(sprintf("Average number of DEGs per method: %.2f", mean(sizes)))
    }

    cons <- consensus(object)
    if (length(cons) == 0L) {
        message("Consensus list is empty or not defined.")
    } else {
        message(sprintf("Number of genes in consensus: %d", length(cons)))
    }

    invisible(NULL)
})
