# Run differential expression analysis with KnowSeq

Executes a two-group differential expression analysis using functions
from the KnowSeq package.

## Usage

``` r
runKnowSeq(
  count,
  groupName,
  numberReplic,
  filterId = "ensembl_gene_id",
  notSapiens = FALSE
)
```

## Arguments

- count:

  A matrix or data.frame of raw counts or abundance values, with genes
  in rows and samples in columns.

- groupName:

  A character vector with exactly two group labels.

- numberReplic:

  A single positive integer indicating the number of replicates per
  group.

- filterId:

  A single character string specifying the annotation filter used by
  [`KnowSeq::getGenesAnnotation()`](https://rdrr.io/pkg/KnowSeq/man/getGenesAnnotation.html).

- notSapiens:

  A single logical value indicating whether a non-human annotation
  dataset should be used.

## Value

A data.frame containing the differential expression results returned by
[`KnowSeq::DEGsExtraction()`](https://rdrr.io/pkg/KnowSeq/man/DEGsExtraction.html).

## Details

This function wraps the main KnowSeq workflow steps:

1.  retrieve gene annotation with
    [`KnowSeq::getGenesAnnotation()`](https://rdrr.io/pkg/KnowSeq/man/getGenesAnnotation.html);

2.  compute expression values with
    [`KnowSeq::calculateGeneExpressionValues()`](https://rdrr.io/pkg/KnowSeq/man/calculateGeneExpressionValues.html);

3.  extract differentially expressed genes with
    [`KnowSeq::DEGsExtraction()`](https://rdrr.io/pkg/KnowSeq/man/DEGsExtraction.html).

The function expects exactly two groups. Row names in `count` must
contain gene identifiers compatible with the selected `filterId`.

Because this workflow depends on external package behavior and may
require internet access, errors raised by internal KnowSeq, cqn, or
mclust calls are rethrown with additional context.

## See also

[`runExpression`](https://costasilvati.github.io/consexpressionR/reference/runExpression.md),
[`expressionDefinition`](https://costasilvati.github.io/consexpressionR/reference/expressionDefinition.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(gse95077)
treats <- c("BM", "JJ")
res <- runKnowSeq(
  count = gse95077,
  groupName = treats,
  numberReplic = 3
)
head(res)
} # }
```
