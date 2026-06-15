# Read table count File and creates a matrix from the given content file.

Read table count File and creates a matrix from the given content file.

## Usage

``` r
readCountFile(tableCountPath = "data/table_count_df.csv", split = ",")
```

## Arguments

- tableCountPath:

  the name of the file which the data are to be read from. Each row of
  the table appears as one line of the file. If it does not contain an
  absolute path, the file name is relative to the current working
  directory, getwd(). This can be a compressed file (see file).
  Alternatively, file can be a readable text-mode connection (which will
  be opened for reading if necessary, and if so closed (and hence
  destroyed) at the end of the function call).

- split:

  the field separator character. Values on each line of the file are
  separated by this character. If sep = "" (the default for read.table)
  the separator is ‘,’, newlines or carriage returns.

## Value

content of file in matrix R format

## Examples

``` r
# Example using a small toy count table
tmp <- tempfile(fileext = ".csv")
write.csv(
  data.frame(
    gene = c("gene1", "gene2"),
    sample1 = c(10, 20),
    sample2 = c(15, 25)
  ),
  tmp,
  row.names = FALSE
)
counts <- readCountFile(tmp)
head(counts)
#>       sample1 sample2
#> gene1      10      15
#> gene2      20      25
```
