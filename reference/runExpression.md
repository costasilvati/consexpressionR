# Make expression analysis of multiple tools and return results by tool.

Make expression analysis of multiple tools and return results by tool.

## Usage

``` r
runExpression(
  numberReplics,
  groupName,
  tableCountPath = NULL,
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
res <- runExpression(
  numberReplics = 2,
  groupName = treats,
  rDataFrameCount = counts,
  controlDeseq2 = "control",
  contrastDeseq2 = "treated",
  printResults = FALSE,
  outDirPath = tempdir()
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
summary(res)
#> Summary of ExpressionResultSet
#> ---------------------------------
#> Methods used: edger, limma, noiseq, ebseq, deseq2, samseq
#> Number of result tables: 6
#> Average number of DEGs per method: 86.67
#> Consensus list is empty or not defined.
```
