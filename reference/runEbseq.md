# Execute EBSeq gene Expression analysis

Execute EBSeq gene Expression analysis

## Usage

``` r
runEbseq(
  countMatrix,
  designExperiment,
  fdr = 0.05,
  ppThreshold = 0.8,
  maxRound = 50,
  methodDeResults = "robust",
  groups = c("")
)
```

## Arguments

- countMatrix:

  either a matrix of raw (read) counts.

- designExperiment:

  replicate and treatment by samples

- fdr:

  parameter used in EBTest function: fdr False Discovery Rate cutt off

- ppThreshold:

  posterior Probability Threshold

- maxRound:

  parameter used in EBTest function: Number of iterations. The default
  value is 50.

- methodDeResults:

  parameter used in GetDEResults function: "robust" or "classic". Using
  the "robust" option, EBSeq is more robust to genes with outliers and
  genes with extremely small variances. Using the "classic" option, the
  results will be more comparable to those obtained by using the
  GetPPMat() function from earlier version (\<= 1.7.0) of EBSeq

- groups:

  text, name of samples or treatment

## Value

EBSeq report in data Frame fromat

## Examples

``` r
data(gse95077)
treats <- c("BM", "JJ")
designExperimentModel <- rep(treats, each = 3)
toolResult <- NULL
toolResult$ebseq <- runEbseq(gse95077,designExperimentModel,groups = treats)
#> Warning: `expect_is()` was deprecated in the 3rd edition.
#> ℹ Use `expect_type()`, `expect_s3_class()`, or `expect_s4_class()` instead
#> Warning: `expect_is()` was deprecated in the 3rd edition.
#> ℹ Use `expect_type()`, `expect_s3_class()`, or `expect_s4_class()` instead
#> Warning: `expect_is()` was deprecated in the 3rd edition.
#> ℹ Use `expect_type()`, `expect_s3_class()`, or `expect_s4_class()` instead
#> Initial number of DE patterns = 2
#> Final number of DE patterns = 2
```
