# Pick power to fit network to a scale-free topology

Pick power to fit network to a scale-free topology

## Usage

``` r
SFT_fit(exp, net_type = "signed", rsquared = 0.8, cor_method = "spearman")
```

## Arguments

- exp:

  A gene expression data frame with genes in row names and samples in
  column names or a \`SummarizedExperiment\` object.

- net_type:

  Network type. One of 'signed', 'signed hybrid' or 'unsigned'. Default
  is signed.

- rsquared:

  R squared cutoff. Default is 0.8.

- cor_method:

  Correlation method. One of "pearson", "biweight" or "spearman".
  Default is "spearman".

## Value

A list containing:

- **power** Numeric, optimal power based on scale-free topology fit

- **plot** A ggplot object displaying main statistics of the SFT fit
  test

## See also

[`pickSoftThreshold`](https://rdrr.io/pkg/WGCNA/man/pickSoftThreshold.html)

## Author

Fabricio Almeida-Silva

## Examples

``` r
data(filt.se)
sft <- SFT_fit(filt.se, cor_method = "pearson")
#> Warning: executing %dopar% sequentially: no parallel backend registered
#>    Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
#> 1      3   0.3340  1.0300        0.18900   145.0     149.0  222.0
#> 2      4   0.2060  0.5420        0.05910   116.0     114.0  195.0
#> 3      5   0.1460  0.3930       -0.00703    96.8      91.3  176.0
#> 4      6   0.0556  0.1930       -0.07170    83.6      74.2  162.0
#> 5      7   0.0209  0.1050       -0.07690    73.7      61.8  151.0
#> 6      8   0.0015 -0.0255       -0.12100    66.0      51.9  142.0
#> 7      9   0.0443 -0.1400       -0.16000    59.8      46.5  134.0
#> 8     10   0.1470 -0.2570       -0.08480    54.7      43.2  127.0
#> 9     11   0.2570 -0.3460        0.04810    50.4      39.7  121.0
#> 10    12   0.3320 -0.3770        0.14800    46.7      36.6  116.0
#> 11    13   0.4520 -0.4400        0.29700    43.5      32.8  111.0
#> 12    14   0.5470 -0.4910        0.41800    40.8      30.1  107.0
#> 13    15   0.5880 -0.5240        0.47300    38.3      27.6  103.0
#> 14    16   0.6870 -0.5470        0.59800    36.1      24.8   99.5
#> 15    17   0.7570 -0.5730        0.68900    34.1      22.9   96.1
#> 16    18   0.8060 -0.5880        0.75100    32.3      21.4   93.0
#> 17    19   0.8350 -0.6170        0.78800    30.7      20.1   90.1
#> 18    20   0.8770 -0.6480        0.84100    29.2      19.3   87.4
```
