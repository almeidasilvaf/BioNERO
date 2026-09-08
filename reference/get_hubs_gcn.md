# Get GCN hubs

Get GCN hubs

## Usage

``` r
get_hubs_gcn(exp, net)
```

## Arguments

- exp:

  A gene expression data frame with genes in row names and samples in
  column names or a \`SummarizedExperiment\` object.

- net:

  List object returned by `exp2gcn`.

## Value

Data frame containing gene IDs, modules and intramodular connectivity of
all hubs.

## See also

[`signedKME`](https://rdrr.io/pkg/WGCNA/man/signedKME.html)

## Author

Fabricio Almeida-Silva

## Examples

``` r
data(filt.se)
gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson")
#> ..connectivity..
#> ..matrix multiplication (system BLAS)..
#> ..normalization..
#> ..done.
hubs <- get_hubs_gcn(filt.se, gcn)
```
