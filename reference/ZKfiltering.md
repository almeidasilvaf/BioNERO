# Filter outlying samples based on the standardized connectivity (Zk) method

Filter outlying samples based on the standardized connectivity (Zk)
method

## Usage

``` r
ZKfiltering(exp, zk = -2, cor_method = "spearman")
```

## Arguments

- exp:

  A gene expression data frame with genes in row names and samples in
  column names or a \`SummarizedExperiment\` object.

- zk:

  Standardized connectivity threshold. Default is -2.

- cor_method:

  Correlation method. One of "pearson", "biweight" or "spearman".
  Default is "spearman".

## Value

Filtered gene expression data frame or \`SummarizedExperiment\` object.

## References

Oldham, M. C., Langfelder, P., & Horvath, S. (2012). Network methods for
describing sample relationships in genomic datasets: application to
Huntington’s disease. BMC systems biology, 6(1), 1-18.

## See also

[`adjacency`](https://rdrr.io/pkg/WGCNA/man/adjacency.html)

## Author

Fabricio Almeida-Silva

## Examples

``` r
data(zma.se)
filt_exp <- ZKfiltering(zma.se)
#> Number of removed samples: 1
```
