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
data(gse95077)
treats <- c("BM", "JJ")
cons_result <- runExpression(numberReplics = 3,
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
expDef_result <- expressionDefinition(cons_result, treats)
deByTool <- listDeByTool(cons_result,rownames(gse95077), expDef_result)

#      DESeq2 limma edgeR
# gene1      1     1     1
# gene2      0     0     1
# gene3      1     0     1
# gene4      0     1     0
# gene5      0     0     0
# gene6      1     1     1
```
