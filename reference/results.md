# Get result tables

Get result tables

## Usage

``` r
results(object)
```

## Arguments

- object:

  An ExpressionResultSet object.

## Value

A named list of data.frames (or NULL elements).

## Examples

``` r
obj <- createExpressionResultSet(
  results = list(edgeR = data.frame(gene = "g1", logFC = 1)),
  methodNames = "edgeR"
)
results(obj)
#> $edgeR
#>   gene logFC
#> 1   g1     1
#> 
```
