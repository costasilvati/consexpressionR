#' @title NCBI GEO access GSE95077 dataset
#'
#' @description Amiloride, an old diuretic drug, is a potential therapeutic agent for multiple myeloma
#' Report ...
#'
#' @format ## `.csv`
#' A data frame with 31 rows and 6 columns:
#' \describe{
#'   \item{design}{Poly A+ RNA from KMS12-BM and JJN3 cells untreated or treated with amiloride or TG003 (0.1 mM, 0.4 mM, and 0.4 mM respectively) for 24 h was isolated and prepared for RNA-seq.}
#'   \item{organism}{Homo Sapiens}
#'   \item{year}{2018}
#'   ...
#' }
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95077>
#' @export
gse95077 <- load("data/gse95077.rda")
