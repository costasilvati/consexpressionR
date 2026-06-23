# listDeByTool: creates a data.frame with genes (rows) and DE tools (columns). Each tool that identifies a gene as DE is marked with 1.

listDeByTool: creates a data.frame with genes (rows) and DE tools
(columns). Each tool that identifies a gene as DE is marked with 1.

## Usage

``` r
listDeByTool(consexpressionList, geneNames, deList)
```

## Arguments

- consexpressionList:

  An `ExpressionResultSet` object with DE results.

- geneNames:

  Character vector. Names of all genes in the experiment.

- deList:

  List. Each element contains a character vector with DE gene names
  (from each tool).

## Value

A data.frame with genes as rows and tools as columns; 1 indicates DE by
that tool.

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
cons_result <- runExpression(numberReplics = 2,
                             groupName = treats,
                             rDataFrameCount = counts,
                             controlDeseq2 = "control",
                             contrastDeseq2 = "treated" )
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
expDef_result <- expressionDefinition(cons_result, treats)
deByTool <- listDeByTool(cons_result,rownames(counts), expDef_result)
```
