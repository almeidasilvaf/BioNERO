# Calculate network preservation between two expression data sets

Calculate network preservation between two expression data sets

## Usage

``` r
module_preservation(
  explist,
  ref_net = NULL,
  test_net = NULL,
  algorithm = "netrep",
  nPerm = 1000,
  nThreads = 1
)
```

## Arguments

- explist:

  List of SummarizedExperiment objects or expression data frames with
  genes (or orthogroups) in row names and samples in column names.

- ref_net:

  Reference network object returned by the function `exp2gcn`.

- test_net:

  Test network object returned by the function `exp2gcn`.

- algorithm:

  Module preservation algorithm to be used. One of 'netrep' (default,
  permutation-based) or WGCNA.

- nPerm:

  Number of permutations. Default: 1000

- nThreads:

  Number of threads to be used for parallel computing. Default: 1

## Value

A list containing the preservation statistics (netrep) or a ggplot
object with preservation statistics. See
[`WGCNA::modulePreservation`](https://rdrr.io/pkg/WGCNA/man/modulePreservation.html)
or
[`NetRep::modulePreservation`](https://rdrr.io/pkg/NetRep/man/modulePreservation.html)
for more info.

## Examples

``` r
# \donttest{
set.seed(1)
data(og.zma.osa)
data(zma.se)
data(osa.se)
og <- og.zma.osa
exp_ortho <- exp_genes2orthogroups(explist, og, summarize = "mean")
#> Error: object 'explist' not found
exp_ortho <- lapply(exp_ortho, function(x) filter_by_variance(x, n=1500))
#> Error: object 'exp_ortho' not found
# Previously calculated SFT powers
powers <- c(13, 15)
gcn_osa <- exp2gcn(exp_ortho$osa, net_type = "signed hybrid",
                   SFTpower = powers[1], cor_method = "pearson")
#> Error: object 'exp_ortho' not found
gcn_zma <- exp2gcn(exp_ortho$zma, net_type = "signed hybrid",
                   SFTpower = powers[2], cor_method = "pearson")
#> Error: object 'exp_ortho' not found
explist <- exp_ortho
#> Error: object 'exp_ortho' not found
ref_net <- gcn_osa
#> Error: object 'gcn_osa' not found
test_net <- gcn_zma
#> Error: object 'gcn_zma' not found
# 10 permutations for demonstration purposes
pres <- module_preservation(explist, ref_net, test_net, nPerm=10)
#> Error: object 'explist' not found
# }
```
