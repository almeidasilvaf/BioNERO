# Perform module stability analysis

Perform module stability analysis

## Usage

``` r
module_stability(exp, net, nRuns = 20)
```

## Arguments

- exp:

  A gene expression data frame with genes in row names and samples in
  column names or a \`SummarizedExperiment\` object.

- net:

  List object returned by `exp2gcn`.

- nRuns:

  Number of times to resample. Default is 20.

## Value

A base plot with the module stability results.

## See also

[`sampledBlockwiseModules`](https://rdrr.io/pkg/WGCNA/man/sampledBlockwiseModules.html)

## Examples

``` r
data(filt.se)
filt <- filt.se[1:100, ] # reducing even further for testing purposes
# The SFT fit was previously calculated and the optimal power was 16
gcn <- exp2gcn(filt, SFTpower = 16, cor_method = "pearson")
#> ..connectivity..
#> ..matrix multiplication (system BLAS)..
#> ..normalization..
#> ..done.
# For simplicity, only 2 runs
module_stability(exp = filt, net = gcn, nRuns = 2)
#>  ...working on run 1 ..
#>  ...working on run 2 ..
#>  ...working on run 3 ..
```
