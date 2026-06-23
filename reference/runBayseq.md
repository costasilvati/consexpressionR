# Run differential expression analysis using baySeq

Performs differential expression analysis on RNA-seq count data using
the baySeq package, which applies empirical Bayesian methods to estimate
posterior likelihoods of differential expression.

## Usage

``` r
runBayseq(
  data,
  groups = NULL,
  sampleInfo = NULL,
  pCutoff = 0.95,
  fdrCutoff = 0.05,
  samplesize = 10000,
  seed = 42
)
```

## Arguments

- data:

  A numeric matrix of raw read counts, with genes as rows and samples as
  columns. Row names should contain gene identifiers.

- groups:

  A named list defining the group structure for baySeq. For a two-group
  comparison, use: `list(NDE = c(1,1,1,1), DE = c(1,1,2,2))`, where
  values indicate group membership for each sample. If `NULL`, groups
  are inferred from `sampleInfo`.

- sampleInfo:

  A character or factor vector indicating group membership for each
  sample (e.g., `c("control","control","treated","treated")`). Used to
  build `groups` automatically when `groups = NULL`.

- pCutoff:

  Numeric. Posterior probability cutoff to call a gene as differentially
  expressed. Default is `0.95`.

- fdrCutoff:

  Numeric. FDR cutoff applied to adjusted p-values for secondary
  filtering. Default is `0.05`.

- samplesize:

  Integer. Number of samples used to estimate priors in baySeq. Default
  is `10000`. Reduce for faster execution during testing.

- seed:

  Integer. Random seed for reproducibility. Default is `42`.

## Value

A data frame with one row per gene, containing:

- gene:

  Gene identifier (from row names of `data`)

- posterior_DE:

  Posterior probability of differential expression

- FDR:

  False discovery rate estimate

- logFC:

  Log2 fold-change between groups (mean group 2 / mean group 1)

- DE_baySeq:

  Logical. `TRUE` if the gene is called as differentially expressed
  based on `pCutoff` and `fdrCutoff`

Returns `NULL` invisibly if baySeq is not installed or if the analysis
fails.

## Details

The function wraps the baySeq workflow for two-group RNA-seq
comparisons. It constructs a `countData` object, estimates prior
distributions via `getPriors.NB`, and computes likelihoods via
`getLikelihoods`. Genes are called as differentially expressed when
their posterior probability of DE exceeds `pCutoff` AND their FDR is
below `fdrCutoff`.

The log2 fold-change is calculated from the group means of the raw count
matrix and is provided for reference only; baySeq itself does not model
fold-change directly.

## References

Hardcastle TJ and Kelly KA (2010). baySeq: Empirical Bayesian methods
for identifying differential expression in sequence count data. BMC
Bioinformatics, 11, 422.
[doi:10.1186/1471-2105-11-422](https://doi.org/10.1186/1471-2105-11-422)

## See also

[`countData`](https://rdrr.io/pkg/baySeq/man/baySeq-classes.html),
[`getPriors.NB`](https://rdrr.io/pkg/baySeq/man/getPriors.html),
[`getLikelihoods`](https://rdrr.io/pkg/baySeq/man/getLikelihoods.html)

## Examples

``` r
if (requireNamespace("baySeq", quietly = TRUE)) {
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
  result <- runBayseq(counts, sampleInfo = groups_info)
  utils::head(result)
}
#> baySeq libsizes
#> baySeq getPriors
#> Finding priors...
#> Warning: The '@replicates' slot is not a factor; converting now.
#> done.
#> Finding posterior likelihoods...
#> Length of priorReps:0
#> Length of priorSubset:100
#> Length of subset:100
#> Length of postRows:100
#> Analysing part 1 of 1
#> Preparing data...
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> .
#> done.
#> Estimating likelihoods...
#> ...done!
#> .
#> done.
#> Available columns in results_raw:
#> [1] "sample1" "sample2" "sample3" "sample4" "likes"   "DE"      "FDR.DE" 
#> [8] "FWER.DE"
#>      sample1 sample2 sample3 sample4     likes  DE    FDR.DE   FWER.DE
#> NA         6       2     163     183 0.7791478 2>1 0.2208522 0.2208522
#> NA.1       0       1      36      39 0.5303252 2>1 0.3452635 0.5867983
#> NA.2       3       3      95     126 0.5212935 2>1 0.3897445 0.7846006
#> NA.3      14      15       4      29 0.5017534 1>2 0.4168700 0.8919226
#> NA.4      32       8      23      29 0.4885126 1>2 0.4357935 0.9472029
#> NA.5       1       2     187     104 0.4878071 2>1 0.4485267 0.9742452
#> Warning: ‘>=’ not meaningful for factors
#>   gene posterior_DE       FDR logFC DE_baySeq
#> 1   NA          2>1 0.2208522    NA     FALSE
#> 2 NA.1          2>1 0.3452635    NA     FALSE
#> 3 NA.2          2>1 0.3897445    NA     FALSE
#> 4 NA.3          1>2 0.4168700    NA     FALSE
#> 5 NA.4          1>2 0.4357935    NA     FALSE
#> 6 NA.5          2>1 0.4485267    NA     FALSE
```
