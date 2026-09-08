# Check scale-free topology fit for a given network

Check scale-free topology fit for a given network

## Usage

``` r
check_SFT(edgelist, net_type = "gcn")
```

## Arguments

- edgelist:

  Edge list as a data frame containing node 1, node 2 and edge weight.

- net_type:

  Type of biological network. One of "gcn", "grn", or "ppi". Default:
  gcn.

## Value

A list with SFT fit statistics and a message indicating if the network
is scale-free.

## Examples

``` r
set.seed(1)
exp <- t(matrix(rnorm(10000), ncol=1000, nrow=200))
rownames(exp) <- paste0("Gene", 1:nrow(exp))
colnames(exp) <- paste0("Sample", 1:ncol(exp))
cormat <- cor(t(exp))
edges <- cormat_to_edgelist(cormat)
edges <- edges[abs(edges$Weight) > 0.10, ]
check_SFT(edges)
#> Could not obtain P-value for the Kolmogorov-Smirnov statistic.
#> $continuous
#> [1] FALSE
#> 
#> $alpha
#> [1] 9.326825
#> 
#> $xmin
#> [1] 219
#> 
#> $logLik
#> [1] -877.4942
#> 
#> $KS.stat
#> [1] 0.1255631
#> 
```
