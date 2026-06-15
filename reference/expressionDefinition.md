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
data(gse95077)
treats <- c("BM", "JJ")
res <- runExpression(numberReplics = 3,
                     groupName = treats,
                     rDataFrameCount = gse95077,
                     controlDeseq2 = "BM",
                     contrastDeseq2 = "JJ" )
#> calcNormFactors has been renamed to normLibSizes
#> Using classic mode.
#> Getting annotation of the Homo Sapiens...
#> Using reference genome 38.
#> Calculating gene expression values...
#> RQ fit ......
#> SQN 
#> Warning: KnowSeq failed and will be skipped: KnowSeq execution failed due to an internal dependency issue involving the 'mclust' package: could not find function "mclustBIC"
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
#> perm= 25
#> perm= 26
#> perm= 27
#> perm= 28
#> perm= 29
#> perm= 30
#> perm= 31
#> perm= 32
#> perm= 33
#> perm= 34
#> perm= 35
#> perm= 36
#> perm= 37
#> perm= 38
#> perm= 39
#> perm= 40
#> perm= 41
#> perm= 42
#> perm= 43
#> perm= 44
#> perm= 45
#> perm= 46
#> perm= 47
#> perm= 48
#> perm= 49
#> perm= 50
#> perm= 51
#> perm= 52
#> perm= 53
#> perm= 54
#> perm= 55
#> perm= 56
#> perm= 57
#> perm= 58
#> perm= 59
#> perm= 60
#> perm= 61
#> perm= 62
#> perm= 63
#> perm= 64
#> perm= 65
#> perm= 66
#> perm= 67
#> perm= 68
#> perm= 69
#> perm= 70
#> perm= 71
#> perm= 72
#> perm= 73
#> perm= 74
#> perm= 75
#> perm= 76
#> perm= 77
#> perm= 78
#> perm= 79
#> perm= 80
#> perm= 81
#> perm= 82
#> perm= 83
#> perm= 84
#> perm= 85
#> perm= 86
#> perm= 87
#> perm= 88
#> perm= 89
#> perm= 90
#> perm= 91
#> perm= 92
#> perm= 93
#> perm= 94
#> perm= 95
#> perm= 96
#> perm= 97
#> perm= 98
#> perm= 99
#> perm= 100
#> Number of thresholds chosen (all possible thresholds) = 8
#> Getting all the cutoffs for the thresholds...
#> Getting number of false positives in the permutation...
#> Differential expression analysis completed.
expDef <- expressionDefinition(res, treats)
```
