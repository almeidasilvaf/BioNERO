# Keep only genes with the highest variances

Keep only genes with the highest variances

## Usage

``` r
filter_by_variance(exp, n = NULL, percentile = NULL)
```

## Arguments

- exp:

  A gene expression data frame with genes in row names and samples in
  column names or a \`SummarizedExperiment\` object.

- n:

  Number of most variable genes (e.g., n=5000 will keep the top 5000
  most variable genes).

- percentile:

  Percentile of most highly variable genes (e.g., percentile=0.1 will
  keep the top 10 percent most variable genes). Values must range from 0
  to 1.

## Value

Expression data frame or \`SummarizedExperiment\` object with the most
variable genes in row names and samples in column names.

## Author

Fabricio Almeida-Silva

## Examples

``` r
data(zma.se)
filt_exp <- filter_by_variance(zma.se, p=0.1)
```
