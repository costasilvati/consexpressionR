# UpSet Plot shows intersections by DE indications

UpSet Plot shows intersections by DE indications

## Usage

``` r
upSetPlotTools(
  df,
  condition = "condition_default",
  pathOut = ".",
  writeData = TRUE
)
```

## Arguments

- df:

  data.frame: output by 'listDeByTool' function.

- condition:

  string: name of analysed experminet (default: "condition_default")

- pathOut:

  string: path to write a pdf file with upset Plot (default: ".")

- writeData:

  boolean: TRUE if want write a PDF file (default: TRUE)

## Value

data.frame with row sum by input df

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
cons_result <- runExpression(2, treats,rDataFrameCount = counts,
                            sepCharacter = ",",experimentName = "test_cons",
                            outDirPath = "." )
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
#> Warning: DESeq2 skipped: controlGroup='' and/or contrastGroup='' not found in groupName.
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
expDef_result <- expressionDefinition(resultTool = cons_result, groups = treats)
deByTool <- listDeByTool(consexpressionList = cons_result,
                        geneNames = rownames(counts),
                        deList = expDef_result)
if (requireNamespace("UpSetR", quietly = TRUE)) {
    upSetPlotData <- upSetPlotTools(df = deByTool,condition = "Control_vs_Treat",
                                pathOut = ".", writeData = FALSE)
}
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue at <https://github.com/hms-dbmi/UpSetR/issues>.
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue at <https://github.com/hms-dbmi/UpSetR/issues>.
#> Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
#> ℹ Please use the `linewidth` argument instead.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue at <https://github.com/hms-dbmi/UpSetR/issues>.
```
