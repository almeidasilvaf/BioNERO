# Remove missing values in a gene expression data frame

Remove missing values in a gene expression data frame

## Usage

``` r
replace_na(exp, replaceby = 0)
```

## Arguments

- exp:

  A gene expression data frame with genes in row names and samples in
  column names or a \`SummarizedExperiment\` object.

- replaceby:

  What to use instead of NAs. One of 0 or 'mean'. Default is 0.

## Value

Gene expression data frame or \`SummarizedExperiment\` object with all
NAs replaced according to the argument 'replaceby'.

## Author

Fabricio Almeida-Silva

## Examples

``` r
data(zma.se)
exp <- replace_na(zma.se)
sum(is.na(exp))
#> [1] 0
```
