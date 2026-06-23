# consexpressionR: Differential Expression Analysis using consensus of multiple methods

Comprehensive workflow for differential expression (DE) analysis using
multiple methods (e.g., DESeq2, edgeR, limma, EBSeq, NOISeq, SAMSeq,
KnowSeq) and consensus strategies to derive robust DE gene lists.

## See also

Useful links:

- <https://costasilvati.github.io/consexpressionR/>

- Report bugs at
  <https://github.com/costasilvati/consexpressionR/issues>

## Author

**Maintainer**: Juliana Costa-Silva <julianacostasilvati@gmail.com>
([ORCID](https://orcid.org/0000-0002-5830-8548)) \[funder\]

Authors:

- Juliana Costa-Silva <julianacostasilvati@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-5830-8548)) \[funder\]

- David Menotti <menottid@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-2430-2030))

- Fabricio M. Lopes <fabricio@utfpr.edu.br>
  ([ORCID](https://orcid.org/0000-0002-8025-6422))

## Examples

``` r
obj <- createExpressionResultSet(
  results = list(edgeR = data.frame(gene = "g1", logFC = 1)),
  methodNames = "edgeR"
)
obj
#> ExpressionResultSet object
#> Methods executed: edgeR
#> Result tables: 1
#> Consensus: not computed or empty.
```
