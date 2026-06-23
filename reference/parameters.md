# Get parameters used in the analysis

Get parameters used in the analysis

## Usage

``` r
parameters(object)

# S4 method for class 'ExpressionResultSet'
parameters(object)
```

## Arguments

- object:

  An ExpressionResultSet object.

## Value

A list.

## Examples

``` r
obj <- createExpressionResultSet(
  results = list(edgeR = data.frame(gene = "g1", logFC = 1)),
  methodNames = "edgeR",
  parameters = list(alpha = 0.01)
)
parameters(obj)
#> $alpha
#> [1] 0.01
#> 
```
