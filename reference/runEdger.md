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
data(gse95077)
treats <- c("BM", "JJ")
numberReplicsModel <- 3
designExperimentModel <- rep(treats, each = numberReplicsModel)
toolResult <- NULL
toolResult$edger <- runEdger(countMatrix = gse95077, numberReplics = numberReplicsModel,
                               desingExperiment = designExperimentModel)
#> calcNormFactors has been renamed to normLibSizes
#> Using classic mode.
```
