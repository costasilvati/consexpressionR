installDependencies <- function (){ #SAMSeq install
    install.packages(c("matrixStats", "GSA", "shiny", "openxlsx", "Rcpp"))
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("impute")
    # install.packages("devtools") if devtools not installed
    #library(devtools)
    #install_github("cran/samr")
    install.packages("samr")
}