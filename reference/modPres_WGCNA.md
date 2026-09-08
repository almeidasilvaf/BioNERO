# Calculate module preservation between two expression data sets using WGCNA's algorithm

Calculate module preservation between two expression data sets using
WGCNA's algorithm

## Usage

``` r
modPres_WGCNA(explist, ref_net, nPerm = 200)
```

## Arguments

- explist:

  List of expression data frames or SummarizedExperiment objects.

- ref_net:

  Reference network object returned by the function `exp2net`.

- nPerm:

  Number of permutations for the module preservation statistics. It must
  be greater than 1. Default: 200.

## Value

A ggplot object with module preservation statistics.

## Examples

``` r
# \donttest{
set.seed(1)
data(og.zma.osa)
data(zma.se)
data(osa.se)
explist <- list(Zma = zma.se, Osa = osa.se)
og <- og.zma.osa
exp_ortho <- exp_genes2orthogroups(explist, og, summarize = "mean")
exp_ortho <- lapply(exp_ortho, function(x) filter_by_variance(x, n=1500))
# Previously calculated power
powers <- c(13, 15)
gcn_osa <- exp2gcn(exp_ortho$Osa, net_type = "signed hybrid",
                   SFTpower = powers[1], cor_method = "pearson")
#> ..connectivity..
#> ..matrix multiplication (system BLAS)..
#> ..normalization..
#> ..done.
explist <- exp_ortho
ref_net <- gcn_osa
# 5 permutations for demonstration purposes
pres_wgcna <- modPres_WGCNA(explist, ref_net, nPerm=5)
# }
```
