# Identify DE genes by method-specific thresholds

Applies specific filtering rules to each expression method result in the
ExpressionResultSet object and returns the DE genes per tool.

## Usage

``` r
expressionDefinition(
  resultTool,
  groups = c(""),
  lfcMinLimma = -2,
  lfcMaxLimma = 2,
  pValueLimma = 0.05,
  FLimma = 0.8,
  lfcMinSamseq = -2,
  lfcMaxSamseq = 2,
  qValueSamseq = 0.05,
  scoreDSamseq = 0.8,
  lfcMinDeseq2 = -2,
  lfcMaxDeseq2 = 2,
  pValueDeseq2 = 0.05,
  lfcMinEdger = -2,
  lfcMaxEdger = 2,
  pValueEdger = 0.05,
  probNoiseq = 0.8,
  lfcMaxKnowseq = 2,
  lfcMinKnowseq = -2,
  pValueKnowseq = 0.05,
  deClassEbseq = "DE",
  ppThresholdEbseq = 0.8,
  printResults = FALSE,
  pathOutput = "."
)
```

## Arguments

- resultTool:

  An `ExpressionResultSet` object.

- groups:

  Character vector of group names or conditions.

- lfcMinLimma:

  Minimum logFC threshold for limma.

- lfcMaxLimma:

  Maximum logFC threshold for limma.

- pValueLimma:

  P-value cutoff for limma.

- FLimma:

  Minimum F-statistic for multi-group limma.

- lfcMinSamseq:

  Minimum fold-change for SAMSeq.

- lfcMaxSamseq:

  Maximum fold-change for SAMSeq.

- qValueSamseq:

  q-value cutoff for SAMSeq.

- scoreDSamseq:

  Score(d) threshold for SAMSeq multi-class.

- lfcMinDeseq2:

  Minimum log2 fold-change for DESeq2.

- lfcMaxDeseq2:

  Maximum log2 fold-change for DESeq2.

- pValueDeseq2:

  P-value cutoff for DESeq2.

- lfcMinEdger:

  Minimum logFC for edgeR.

- lfcMaxEdger:

  Maximum logFC for edgeR.

- pValueEdger:

  P-value cutoff for edgeR.

- probNoiseq:

  Probability threshold for NOISeq.

- lfcMaxKnowseq:

  Maximum logFC for KnowSeq.

- lfcMinKnowseq:

  Minimum logFC for KnowSeq.

- pValueKnowseq:

  P-value cutoff for KnowSeq.

- deClassEbseq:

  DE class to filter for EBSeq (e.g. "DE" or "EE").

- ppThresholdEbseq:

  Posterior probability threshold for EBSeq.

- printResults:

  Logical. Whether to write results to disk.

- pathOutput:

  Character. Directory path to save result tables.

## Value

A list of data.frames with DE genes per method.

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
groups_info <- c("control", "control", "treated", "treated")
treats <- c("control", "treated")
res <- runExpression(
  numberReplics = 2,
  groupName = treats,
  rDataFrameCount = counts,
  controlDeseq2 = "control",
  contrastDeseq2 = "treated"
)
#> calcNormFactors has been renamed to normLibSizes
#> Using classic mode.
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
#> Warning: ‘>=’ not meaningful for factors
#> calcNormFactors has been renamed to normLibSizes
#> Removing intercept from test coefficients
#> [1] "Computing (M,D) values..."
#> [1] "Computing probability of differential expression..."
#> Warning: `expect_is()` was deprecated in the 3rd edition.
#> ℹ Use `expect_type()`, `expect_s3_class()`, or `expect_s4_class()` instead
#> Warning: `expect_is()` was deprecated in the 3rd edition.
#> ℹ Use `expect_type()`, `expect_s3_class()`, or `expect_s4_class()` instead
#> Warning: `expect_is()` was deprecated in the 3rd edition.
#> ℹ Use `expect_type()`, `expect_s3_class()`, or `expect_s4_class()` instead
#> Initial number of DE patterns = 2
#> Final number of DE patterns = 2
#> converting counts to integer mode
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
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
#> Number of thresholds chosen (all possible thresholds) = 11
#> Getting all the cutoffs for the thresholds...
#> Getting number of false positives in the permutation...
#> Differential expression analysis completed.
expDef <- expressionDefinition(res, treats)
```
