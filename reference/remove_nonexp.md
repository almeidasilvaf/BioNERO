# Remove genes that are not expressed based on a user-defined threshold

Remove genes that are not expressed based on a user-defined threshold

## Usage

``` r
remove_nonexp(
  exp,
  method = "median",
  min_exp = 1,
  min_percentage_samples = 0.25
)
```

## Arguments

- exp:

  A gene expression data frame with genes in row names and samples in
  column names or a \`SummarizedExperiment\` object.

- method:

  Criterion to filter non-expressed genes out. One of "mean", "median",
  "percentage", or "allsamples". Default is "median".

- min_exp:

  If method is 'mean', 'median', or 'allsamples', the minimum value for
  a gene to be considered expressed. If method is 'percentage', the
  minimum value each gene must have in at least n percent of samples to
  be considered expressed.

- min_percentage_samples:

  In case the user chooses 'percentage' as method, expressed genes must
  have expression \>= min_exp in at least this percentage. Values must
  range from 0 to 1.

## Value

Filtered gene expression data frame or \`SummarizedExperiment\` object.

## See also

[`rowMedians`](https://rdrr.io/pkg/matrixStats/man/rowMedians.html)
[`goodSamplesGenes`](https://rdrr.io/pkg/WGCNA/man/goodSamplesGenes.html)

## Author

Fabricio Almeida-Silva

## Examples

``` r
data(zma.se)
filt_exp <- remove_nonexp(zma.se, min_exp = 5)
```
