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
data(gse95077)
treats <- c("BM", "JJ")
toolResult <- NULL
toolResult$limma <- runLimma(gse95077, 3, rep(treats, each = 3))
#> calcNormFactors has been renamed to normLibSizes
#> Removing intercept from test coefficients
```
