# Write a table or data.frame in defined output file

Write a table or data.frame in defined output file

## Usage

``` r
writeResults(
  data,
  toolName = "toolDE_x",
  sepCharacter = "\t",
  pathOutput = "."
)
```

## Arguments

- data:

  table or data.frame dataset

- toolName:

  text with file name

- sepCharacter:

  character to separate data columns

- pathOutput:

  path to write output, needs to be a directory (default:".")

## Value

The path of the file created (invisible)

## Examples

``` r
df <- data.frame(
  col1 = c("treat1", "treat2", "treat3"),
  col2 = c(1, 2, 3)
)
out <- writeResults(data = df, toolName = "test", pathOutput = tempdir())
file.exists(out)
#> [1] TRUE
```
