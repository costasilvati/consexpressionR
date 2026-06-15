# Execute NOISeq gene Expression analysis

Execute NOISeq gene Expression analysis

## Usage

``` r
runNoiSeq(
  countMatrix,
  designExperiment,
  groups = c(""),
  normParm = "rpkm",
  kParam = 0.5,
  factorParam = c("Tissue"),
  lcParam = 0,
  replicatesParam = "technical",
  condExp = c("")
)
```

## Arguments

- countMatrix:

  either a matrix of raw (read) counts.

- designExperiment:

  replicate and treatment by samples.

- groups:

  text list separated by comma, name of samples or treatment.

- normParm:

  Normalization method t can be one of "rpkm" (default), "uqua" (upper
  quartile), "tmm" (trimmed mean of M) or "n" (no normalization).

- kParam:

  Counts equal to 0 are replaced by k. By default, k = 0.5

- factorParam:

  A string indicating the name of factor whose levels are the conditions
  to be compared.

- lcParam:

  Length correction is done by dividing expression by length^lc. By
  default, lc = 0

- replicatesParam:

  In this argument, the type of replicates to be used is defined:
  "technical", "biological" or "no" replicates. By default, "technical"
  replicates option is chosen.

- condExp:

  A vector containing the two conditions to be compared by the
  differential expression algorithm (needed when the factor contains
  more than 2 different conditions).

## Value

A data.frame with output of noiseqOut@results

## Examples

``` r
data(gse95077)
treats <- c("BM", "JJ")
toolResult <- NULL
toolResult$noiseq <- runNoiSeq(countMatrix = gse95077, designExperiment = rep(treats, each = 3))
#> [1] "Computing (M,D) values..."
#> [1] "Computing probability of differential expression..."
```
