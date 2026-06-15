# An S4 class to store differential expression results

An S4 class to store differential expression results

## Slots

- `results`:

  A named list of data.frames, each corresponding to a method (e.g.,
  edgeR, DESeq2, etc.)

- `methodNames`:

  A character vector of methods used

- `parameters`:

  A list of parameters used in the analysis

- `consensus`:

  A list containing consensus DEG results

## Examples

``` r
obj <- createExpressionResultSet(
  results = list(edgeR = data.frame(gene = "g1", logFC = 1)),
  methodNames = "edgeR",
  parameters = list(),
  consensus = list()
)
obj
#> ExpressionResultSet object
#> Methods executed: edgeR
#> Result tables: 1
#> Consensus: not computed or empty.
```
