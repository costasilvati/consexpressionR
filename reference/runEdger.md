# Execute edgeR expression Analisys

Execute edgeR expression Analisys

## Usage

``` r
runEdger(countMatrix, numberReplics, desingExperiment, methNorm = "TMM")
```

## Arguments

- countMatrix:

  either a matrix of raw (read) counts.

- numberReplics:

  number of replicate (technical or biologcal) integer

- desingExperiment:

  replicate and treatment by samples

- methNorm:

  normalization method to be used in edgeR::calcNormFactors(), default:
  "TMM"

## Value

edgeR report in data Frame

## Examples

``` r
set.seed(42)
counts <- matrix(
  as.integer(c(
    rnbinom(200, mu = 10, size = 1),
    rnbinom(200, mu = 100, size = 1)
  )),
  nrow = 100,
  dimnames = list(
    paste0("gene", seq_len(100)),
    paste0("sample", seq_len(4))
  )
)
treats <- c("control", "treated")
numberReplicsModel <- 2
designExperimentModel <- rep(treats, each = numberReplicsModel)
toolResult <- NULL
toolResult$edger <- runEdger(
  countMatrix = counts, numberReplics = numberReplicsModel,
  desingExperiment = designExperimentModel
)
#> calcNormFactors has been renamed to normLibSizes
#> Using classic mode.
```
