#' @keywords internal
#' @noRd
.check_package <- function(pkg, repo = "CRAN") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install_cmd <- if (repo == "Bioconductor") {
      sprintf("BiocManager::install('%s')", pkg)
    } else {
      sprintf("install.packages('%s')", pkg)
    }
    stop(sprintf(
      "Package '%s' is necessary to execute this analysis, but is not installed.\n  Install using R command: %s",
      pkg, install_cmd
    ), call. = FALSE)
  }
}

#' @keywords internal
#' @noRd
.progress <- function(detail, progress_obj = NULL) {
  if (!is.null(progress_obj)) {
    progress_obj$set(detail = detail)
  }
}

#' @keywords internal
#' @noRd
.is_count_like <- function(x) {
  x <- as.matrix(x)
  if (anyNA(x) || any(x < 0)) {
    return(FALSE)
  }
  all(abs(x - round(x)) < .Machine$double.eps^0.5)
}
