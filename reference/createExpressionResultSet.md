# Constructor for ExpressionResultSet

Constructor for ExpressionResultSet

## Usage

``` r
createExpressionResultSet(
  results,
  methodNames,
  parameters = list(),
  consensus = list()
)
```

## Arguments

- results:

  Named list of data.frames (or NULL elements).

- methodNames:

  Character vector with method names.

- parameters:

  List of analysis parameters.

- consensus:

  List with consensus results (default: empty list).

## Value

An ExpressionResultSet object.

## Examples

``` r
res_list <- list(
  edgeR = data.frame(gene = c("g1","g2"), logFC = c(1.2, -0.4)),
  DESeq2 = NULL
)

obj <- createExpressionResultSet(
  results = res_list,
  methodNames = c("edgeR", "DESeq2"),
  parameters = list(alpha = 0.05),
  consensus = list(g1 = TRUE)
)

obj
#> ExpressionResultSet object
#> Methods executed: edgeR, DESeq2
#> Result tables: 2
#> Consensus: available (1 gene(s))
summary(obj)
#>              Length               Class                Mode 
#>                   1 ExpressionResultSet                  S4 
```
