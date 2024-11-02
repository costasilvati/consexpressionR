---
title: "consexpressionR"
output: default
vignette: >
  %\VignetteIndexEntry{consexpressionR-vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# consexpressionR

R package consexpressionR is an R package for differential expression analysis of count or quantification data. With just one count or quantification file, the user can obtain a list of genes indicated as differentially expressed by 5 or more methods. With a user-friendly graphical interface, it can be used in a web browser, where the user can parameterize all 7 expression analysis methodologies.

The differential expression analysis is performed with the methods edgeR, DeSeq2, SAMSeq, KnowSeq, EBSeq, limma, and NOISeq. The results of each tool can be written in .csv files. Visualization is carried out through a heat map Plot, containing the genes that were considered as differentially expressed by 5 or more methods. consexpressionR is an R package for differential expression analysis of count or quantification data.

With just one count or quantification file, the user can obtain a list of genes indicated as differentially expressed by 5 or more methods. With a user-friendly graphical interface, it can be used in a web browser, where the user can parameterize all 7 expression analysis methodologies. The differential expression analysis is performed with the methods edgeR, DeSeq2, SAMSeq, KnowSeq, EBSeq, limma, and NOISeq. The results of each tool can be written in .csv files. Visualization is carried out through a radarPlot, containing the genes thatwere considered as differentially expressed by 5 or more methods.

You can see pkgdown in action at <https://costasilvati.shinyapps.io/consexpressionRapp/>.

Learn more in `vignette("consexpressionR")`.

## Installation

::: consexpressionR-release
``` r
# Install released version from github
library(devtools)
install_github("costasilvati/consexpressionR")
library(consexpressionR)
# shinyApp running in localhost
consexpressionR()
```
:::

## Usage - shiny app

To run shiny aplication and view consexpressionR in your browser:

``` r
library(consexpressionR)
# shinyApp running in localhost
consexpressionR()
```
