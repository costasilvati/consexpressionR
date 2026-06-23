# Introduction to cosexpressionR

## Basics

### Install `consexpressionR`

`R` is an open-source statistical environment which can be easily
modified to enhance its functionality via packages.
*[consexpressionR](https://bioconductor.org/packages/3.23/consexpressionR)*
is a `R` package available via the
[Bioconductor](http://bioconductor.org) repository for packages. `R` can
be installed on any operating system from
[CRAN](https://cran.r-project.org/) after which you can install
*[consexpressionR](https://bioconductor.org/packages/3.23/consexpressionR)*
by using the following commands in your `R` session:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("consexpressionR")

## Check that you have a valid Bioconductor installation
BiocManager::valid()
```

### Required knowledge

*[consexpressionR](https://bioconductor.org/packages/3.23/consexpressionR)*
is based on many other packages and in particular in those that have
implemented the infrastructure needed for dealing with RNA-seq data
(EDIT!). That is, packages like
*[SummarizedExperiment](https://bioconductor.org/packages/3.23/SummarizedExperiment)*
(EDIT!).

If you are asking yourself the question “Where do I start using
Bioconductor?” you might be interested in [this blog
post](http://lcolladotor.github.io/2014/10/16/startbioc/#.VkOKbq6rRuU).

### Asking for help

As package developers, we try to explain clearly how to use our packages
and in which order to use the functions. But `R` and `Bioconductor` have
a steep learning curve so it is critical to learn where to ask for help.
The blog post quoted above mentions some but we would like to highlight
the [Bioconductor support site](https://support.bioconductor.org/) as
the main resource for getting help: remember to use the
`consexpressionR` tag and check [the older
posts](https://support.bioconductor.org/tag/consexpressionR/). Other
alternatives are available such as creating GitHub issues and tweeting.
However, please note that if you want to receive help you should adhere
to the [posting
guidelines](http://www.bioconductor.org/help/support/posting-guide/). It
is particularly critical that you provide a small reproducible example
and your session information so package developers can track down the
source of the error.

### Citing `consexpressionR`

We hope that
*[consexpressionR](https://bioconductor.org/packages/3.23/consexpressionR)*
will be useful for your research. Please use the following information
to cite the package and the overall approach. Thank you!

``` r

## Citation info
citation("consexpressionR")
#> To cite package 'consexpressionR' in publications use:
#> 
#>   Costa-Silva, J., Menotti, D., Lopes, M. F. consexpressionR: an R
#>   package for consensus differential gene expression analysis arxiv
#>   (2025)
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {consexpressionR: an R package for consensus differential gene expression analysis},
#>     author = {Juliana Costa-Silva and David Menotti and Fabricio M. Lopes},
#>     year = {2025},
#>     journal = {arxiv},
#>     doi = {10.48550/arXiv.2503.21546},
#>   }
```

## Quick start to using `consexpressionR`

``` r

library("consexpressionR")
```

Edit this as you see fit =)

Here is an example of you can cite your package inside the vignette:

- *[consexpressionR](https://bioconductor.org/packages/3.23/consexpressionR)*
  (Costa-Silva, Menotti, and Lopes, 2025)

## Reproducibility

The
*[consexpressionR](https://bioconductor.org/packages/3.23/consexpressionR)*
package (Costa-Silva, Menotti, and Lopes, 2025) was made possible thanks
to:

- R (R Core Team, 2026)
- *[BiocStyle](https://bioconductor.org/packages/3.23/BiocStyle)* (Oleś,
  2026)
- *[knitr](https://CRAN.R-project.org/package=knitr)* (Xie, 2025)
- *[RefManageR](https://CRAN.R-project.org/package=RefManageR)* (McLean,
  2017)
- *[rmarkdown](https://CRAN.R-project.org/package=rmarkdown)* (Allaire,
  Xie, Dervieux, McPherson, Luraschi, Ushey, Atkins, Wickham, Cheng,
  Chang, and Iannone, 2026)
- *[sessioninfo](https://CRAN.R-project.org/package=sessioninfo)*
  (Wickham, Chang, Flight, Müller, and Hester, 2026)
- *[testthat](https://CRAN.R-project.org/package=testthat)* (Wickham,
  2011)

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.23/biocthis)*.

`R` session information.

    #> ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
    #>  setting  value
    #>  version  R version 4.6.0 (2026-04-24)
    #>  os       Ubuntu 24.04.4 LTS
    #>  system   x86_64, linux-gnu
    #>  ui       X11
    #>  language en
    #>  collate  C.UTF-8
    #>  ctype    C.UTF-8
    #>  tz       UTC
    #>  date     2026-06-23
    #>  pandoc   3.8.3 @ /opt/hostedtoolcache/pandoc/3.8.3/x64/ (via rmarkdown)
    #>  quarto   NA
    #> 
    #> ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
    #>  package         * version   date (UTC) lib source
    #>  backports         1.5.1     2026-04-03 [1] RSPM
    #>  bibtex            0.5.2     2026-02-03 [1] RSPM
    #>  BiocManager       1.30.27   2025-11-14 [1] RSPM
    #>  BiocStyle       * 2.40.0    2026-04-28 [1] Bioconduc~
    #>  bookdown          0.47      2026-06-16 [1] RSPM
    #>  bslib             0.11.0    2026-05-16 [1] RSPM
    #>  cachem            1.1.0     2024-05-16 [1] RSPM
    #>  cli               3.6.6     2026-04-09 [1] RSPM
    #>  consexpressionR * 0.99.0    2026-06-23 [1] local
    #>  desc              1.4.3     2023-12-10 [1] RSPM
    #>  digest            0.6.39    2025-11-19 [1] RSPM
    #>  evaluate          1.0.5     2025-08-27 [1] RSPM
    #>  fastmap           1.2.0     2024-05-15 [1] RSPM
    #>  fs                2.1.0     2026-04-18 [1] RSPM
    #>  generics          0.1.4     2025-05-09 [1] RSPM
    #>  glue              1.8.1     2026-04-17 [1] RSPM
    #>  htmltools         0.5.9     2025-12-04 [1] RSPM
    #>  htmlwidgets       1.6.4     2023-12-06 [1] RSPM
    #>  httr              1.4.8     2026-02-13 [1] RSPM
    #>  jquerylib         0.1.4     2021-04-26 [1] RSPM
    #>  jsonlite          2.0.0     2025-03-27 [1] RSPM
    #>  knitr             1.51      2025-12-20 [1] RSPM
    #>  lifecycle         1.0.5     2026-01-08 [1] RSPM
    #>  lubridate         1.9.5     2026-02-04 [1] RSPM
    #>  magrittr          2.0.5     2026-04-04 [1] RSPM
    #>  otel              0.2.0     2025-08-29 [1] RSPM
    #>  pkgdown           2.2.0     2025-11-06 [1] RSPM
    #>  plyr              1.8.9     2023-10-02 [1] RSPM
    #>  R6                2.6.1     2025-02-15 [1] RSPM
    #>  ragg              1.5.2     2026-03-23 [1] RSPM
    #>  Rcpp              1.1.1-1.1 2026-04-24 [1] RSPM
    #>  RefManageR      * 1.4.0     2022-09-30 [1] RSPM
    #>  rlang             1.2.0     2026-04-06 [1] RSPM
    #>  rmarkdown         2.31      2026-03-26 [1] RSPM
    #>  sass              0.4.10    2025-04-11 [1] RSPM
    #>  sessioninfo     * 1.2.4     2026-06-04 [1] RSPM
    #>  stringi           1.8.7     2025-03-27 [1] RSPM
    #>  stringr           1.6.0     2025-11-04 [1] RSPM
    #>  systemfonts       1.3.2     2026-03-05 [1] RSPM
    #>  textshaping       1.0.5     2026-03-06 [1] RSPM
    #>  timechange        0.4.0     2026-01-29 [1] RSPM
    #>  xfun              0.59      2026-06-19 [1] RSPM
    #>  xml2              1.6.0     2026-06-22 [1] RSPM
    #>  yaml              2.3.12    2025-12-10 [1] RSPM
    #> 
    #>  [1] /home/runner/work/_temp/Library
    #>  [2] /opt/R/4.6.0/lib/R/site-library
    #>  [3] /opt/R/4.6.0/lib/R/library
    #>  * ── Packages attached to the search path.
    #> 
    #> ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

## Bibliography

This vignette was generated using
*[BiocStyle](https://bioconductor.org/packages/3.23/BiocStyle)* (Oleś,
2026) with *[knitr](https://CRAN.R-project.org/package=knitr)* (Xie,
2025) and *[rmarkdown](https://CRAN.R-project.org/package=rmarkdown)*
(Allaire, Xie, Dervieux et al., 2026) running behind the scenes.

Citations made with
*[RefManageR](https://CRAN.R-project.org/package=RefManageR)* (McLean,
2017).

[\[1\]](#cite-allaire2026rmarkdown) J. Allaire, Y. Xie, C. Dervieux, et
al. *rmarkdown: Dynamic Documents for R*. R package version 2.31. 2026.
URL: <https://github.com/rstudio/rmarkdown>.

[\[2\]](#cite-costasilva2025consexpressionr) J. Costa-Silva, D. Menotti,
and F. M. Lopes. “consexpressionR: an R package for consensus
differential gene expression analysis”. In: *arxiv* (2025). DOI:
[10.48550/arXiv.2503.21546](https://doi.org/10.48550/arXiv.2503.21546).

[\[3\]](#cite-mclean2017refmanager) M. W. McLean. “RefManageR: Import
and Manage BibTeX and BibLaTeX References in R”. In: *The Journal of
Open Source Software* (2017). DOI:
[10.21105/joss.00338](https://doi.org/10.21105/joss.00338).

[\[4\]](#cite-ole2026biocstyle) A. Oleś. *BiocStyle: Standard styles for
vignettes and other Bioconductor documents*. R package version 2.40.0.
2026. DOI:
[10.18129/B9.bioc.BiocStyle](https://doi.org/10.18129/B9.bioc.BiocStyle).
URL: <https://bioconductor.org/packages/BiocStyle>.

[\[5\]](#cite-2026language) R Core Team. *R: A Language and Environment
for Statistical Computing*. R Foundation for Statistical Computing (ROR:
\<<https://ror.org/05qewa988%3E>;). Vienna, Austria, 2026. DOI:
[10.32614/R.manuals](https://doi.org/10.32614/R.manuals). URL:
<https://www.R-project.org/>.

[\[6\]](#cite-wickham2011testthat) H. Wickham. “testthat: Get Started
with Testing”. In: *The R Journal* 3 (2011), pp. 5–10. URL:
<https://journal.r-project.org/articles/RJ-2011-002/>.

[\[7\]](#cite-wickham2026sessioninfo) H. Wickham, W. Chang, R. Flight,
et al. *sessioninfo: R Session Information*. R package version 1.2.4.
2026. URL: <https://sessioninfo.r-lib.org>.

[\[8\]](#cite-xie2025knitr) Y. Xie. *knitr: A General-Purpose Package
for Dynamic Report Generation in R*. R package version 1.51. 2025. URL:
<https://yihui.org/knitr/>.
