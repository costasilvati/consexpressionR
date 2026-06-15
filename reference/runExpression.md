# Make expression analysis of multiple tools and return results by tool.

Make expression analysis of multiple tools and return results by tool.

## Usage

``` r
runExpression(
  numberReplics,
  groupName,
  tableCountPath = "data/gse95077.csv",
  sepCharacter = ",",
  rDataFrameCount = NULL,
  experimentName = "genericExperiment",
  outDirPath = tempdir(),
  printResults = FALSE,
  fitTypeDeseq2 = "local",
  controlDeseq2 = "",
  contrastDeseq2 = "",
  methodNormLimma = "TMM",
  methodAdjPvalueLimma = "BH",
  numberTopTableLimma = 1e+06,
  filterIdKnowseq = "ensembl_gene_id",
  notSapiensKnowseq = FALSE,
  methodNormEdgeR = "TMM",
  normNoiseq = "rpkm",
  kNoiseq = 0.5,
  factorNoiseq = "Tissue",
  lcNoiseq = 0,
  replicatesNoiseq = "technical",
  condExpNoiseq = c(""),
  respTypeSamseq = "Two class unpaired",
  npermSamseq = 100,
  fdrEbseq = 0.05,
  maxRoundEbseq = 50,
  methodDeResultsEbseq = "robust",
  deNovoAanalysis = FALSE,
  progressShiny = NULL
)
```

## Arguments

- numberReplics:

  number of replicates (technical or biological) per group

- groupName:

  character vector, group/treatment names (e.g., c("BM","JJ"))

- tableCountPath:

  path to csv file that contains count/abundance data (local)

- sepCharacter:

  separator used to split csv data (comma or tab)

- rDataFrameCount:

  data.frame/matrix with table count; rownames are genes and colnames
  are samples

- experimentName:

  experiment name

- outDirPath:

  output directory (created if needed)

- printResults:

  if TRUE, writes per-tool outputs to outDirPath

- fitTypeDeseq2:

  DESeq2 fitType

- controlDeseq2:

  control group label for DESeq2

- contrastDeseq2:

  contrast group label for DESeq2

- methodNormLimma:

  normalization method used by limma/edgeR

- methodAdjPvalueLimma:

  p-value adjustment method for limma

- numberTopTableLimma:

  max number of genes for limma topTable

- filterIdKnowseq:

  KnowSeq parameter

- notSapiensKnowseq:

  KnowSeq parameter

- methodNormEdgeR:

  normalization method for edgeR

- normNoiseq:

  NOISeq parameter

- kNoiseq:

  NOISeq parameter

- factorNoiseq:

  NOISeq parameter

- lcNoiseq:

  NOISeq parameter

- replicatesNoiseq:

  NOISeq parameter

- condExpNoiseq:

  NOISeq parameter

- respTypeSamseq:

  SAMSeq parameter

- npermSamseq:

  SAMSeq parameter

- fdrEbseq:

  EBSeq parameter

- maxRoundEbseq:

  EBSeq parameter

- methodDeResultsEbseq:

  EBSeq parameter

- deNovoAanalysis:

  if TRUE, skips KnowSeq (requires annotation)

- progressShiny:

  optional shiny progress callback

## Value

An `ExpressionResultSet` object containing results per method.

## Examples

``` r
data(gse95077)
treats <- c("BM", "JJ")
res <- runExpression(
  numberReplics = 3,
  groupName = treats,
  rDataFrameCount = gse95077,
  controlDeseq2 = "BM",
  contrastDeseq2 = "JJ",
  printResults = FALSE,
  outDirPath = tempdir()
)
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
#> Number of thresholds chosen (all possible thresholds) = 7
#> Getting all the cutoffs for the thresholds...
#> Getting number of false positives in the permutation...
#> Differential expression analysis completed.
summary(res)
#>              Length               Class                Mode 
#>                   1 ExpressionResultSet                  S4 
```
