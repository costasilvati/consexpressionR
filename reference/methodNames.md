# Get method names used in the analysis

Get method names used in the analysis

## Usage

``` r
methodNames(object)

# S4 method for class 'ExpressionResultSet'
methodNames(object)
```

## Arguments

- object:

  An ExpressionResultSet object.

## Value

A character vector.

## Examples

``` r
obj <- createExpressionResultSet(
  results = list(edgeR = data.frame(gene = "g1", logFC = 1)),
  methodNames = c("edgeR")
)
methodNames(obj)
#> [1] "edgeR"
```
