# Execute edgeR expression Analisys

Execute edgeR expression Analisys

## Usage

``` r
runLimma(
  countMatrix,
  numberReplics,
  designExperiment,
  methodNorm = "TMM",
  methodAdjPvalue = "BH",
  numberTopTable = 1e+06
)
```

## Arguments

- countMatrix:

  either a matrix of raw (read) counts.

- numberReplics:

  number of replicate (technical or biologcal) integer

- designExperiment:

  replicate and treatment by samples

- methodNorm:

  normalization method to be used

- methodAdjPvalue:

  correction method, a character string. Can be abbreviated.

- numberTopTable:

  maximum number of genes to list

## Value

limma report in data Frame fromat

## Examples

``` r
set.seed(42)
counts <- matrix(
  as.integer(c(
    rnbinom(200, mu = 10,  size = 1), 
    rnbinom(200, mu = 100, size = 1) 
  )),
  nrow = 100,
  dimnames = list(paste0("gene", seq_len(100)),
                    paste0("sample", seq_len(4)))
)
groups_info <- c("control", "control", "treated", "treated")
treats <- c("control", "treated")
toolResult <- NULL
toolResult$limma <- runLimma(counts, 2, rep(treats, each = 2))
#> calcNormFactors has been renamed to normLibSizes
#> Removing intercept from test coefficients
```
