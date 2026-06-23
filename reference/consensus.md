# Get consensus results

Get consensus results

## Usage

``` r
consensus(object)

# S4 method for class 'ExpressionResultSet'
consensus(object)
```

## Arguments

- object:

  An ExpressionResultSet object.

## Value

A list.

## Examples

``` r
obj <- createExpressionResultSet(
  results = list(edgeR = data.frame(gene = c("g1","g2"), logFC = c(1,-1))),
  methodNames = "edgeR",
  consensus = list(g1 = TRUE)
)
consensus(obj)
#> $g1
#> [1] TRUE
#> 
```
