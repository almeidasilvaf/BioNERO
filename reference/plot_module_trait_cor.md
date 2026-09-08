# Plot a heatmap of module-trait correlations

Plot a heatmap of module-trait correlations

## Usage

``` r
plot_module_trait_cor(corandp, palette = "RdYlBu", transpose = FALSE)
```

## Arguments

- corandp:

  A data frame of module-trait correlations as returned by
  [`module_trait_cor()`](module_trait_cor.md).

- palette:

  Character indicating which RColorBrewer palette to use. Default:
  'RdYlBu'.

- transpose:

  Logical indicating whether to transpose the heatmap or not.

## Value

A \`Heatmap\` object created by
[`ComplexHeatmap::pheatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/pheatmap.html).

## Details

Significance levels: 1 asterisk: significant at alpha = 0.05. 2
asterisks: significant at alpha = 0.01. 3 asterisks: significant at
alpha = 0.001. no asterisk: not significant.

## Examples

``` r
data(filt.se)
gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson")
#> ..connectivity..
#> ..matrix multiplication (system BLAS)..
#> ..normalization..
#> ..done.
corandp <- module_trait_cor(filt.se, MEs = gcn$MEs)
plot_module_trait_cor(corandp)
```
