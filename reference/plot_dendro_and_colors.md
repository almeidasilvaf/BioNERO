# Plot dendrogram of genes and modules

Plot dendrogram of genes and modules

## Usage

``` r
plot_dendro_and_colors(gcn)
```

## Arguments

- gcn:

  List object returned by `exp2gcn`.

## Value

A base plot with the gene dendrogram and modules.

## Examples

``` r
data(filt.se)
gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson")
#> ..connectivity..
#> ..matrix multiplication (system BLAS)..
#> ..normalization..
#> ..done.
plot_dendro_and_colors(gcn)
```
