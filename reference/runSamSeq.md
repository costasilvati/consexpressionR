# Execute SAMSeq gene Expression analysis to count data

Execute SAMSeq gene Expression analysis to count data

## Usage

``` r
runSamSeq(
  countMatrix,
  designExperiment,
  respType = "Two class unpaired",
  numberPermutations = 100
)
```

## Arguments

- countMatrix:

  either a matrix of raw (read) counts.

- designExperiment:

  replicate and treatment by samples

- respType:

  Problem type: "Quantitative" for a continuous parameter; "Two class
  unpaired" for two classes with unpaired observations; "Survival" for
  censored survival outcome; "Multiclass": more than 2 groups; "Two
  class paired" for two classes with paired observations.

- numberPermutations:

  Number of permutations used to estimate false discovery rates

## Value

SAMSeq report in data Frame

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
if (requireNamespace("samr", quietly = TRUE)) {
  toolResult$samseq <- runSamSeq(countMatrix = counts,
                               designExperiment = rep(treats, each = 2))
}
#> Estimating sequencing depths...
#> Resampling to get new data matrices...
#> perm= 1
#> perm= 2
#> perm= 3
#> perm= 4
#> perm= 5
#> perm= 6
#> perm= 7
#> perm= 8
#> perm= 9
#> perm= 10
#> perm= 11
#> perm= 12
#> perm= 13
#> perm= 14
#> perm= 15
#> perm= 16
#> perm= 17
#> perm= 18
#> perm= 19
#> perm= 20
#> perm= 21
#> perm= 22
#> perm= 23
#> perm= 24
#> Number of thresholds chosen (all possible thresholds) = 12
#> Getting all the cutoffs for the thresholds...
#> Getting number of false positives in the permutation...
```
