# Pick power to fit networks to scale-free topology

Pick power to fit networks to scale-free topology

## Usage

``` r
consensus_SFT_fit(
  exp_list,
  setLabels = NULL,
  metadata = NULL,
  cor_method = "spearman",
  net_type = "signed hybrid",
  rsquared = 0.8
)
```

## Arguments

- exp_list:

  A list of expression data frames or SummarizedExperiment objects. If
  input is a list of data frames, row names must correspond to gene IDs
  and column names to samples. The list can be created with
  `list(exp1, exp2, ..., expn)`.

- setLabels:

  Character vector containing labels for each expression set.

- metadata:

  A data frame containing sample names in row names and sample
  annotation in the first column. Ignored if \`exp_list\` is a list of
  \`SummarizedExperiment\` objects, since the function will extract
  colData.

- cor_method:

  Correlation method used for network reconstruction. One of "spearman"
  (default), "biweight", or "pearson".

- net_type:

  Network type. One of "signed hybrid" (default), "signed" or
  "unsigned".

- rsquared:

  Minimum R squared to consider the network similar to a scale-free
  topology. Default is 0.8.

## Value

A list of 2 elements:

- power:

  Numeric vector of optimal beta powers to fit networks to SFT

- plot:

  A ggplot object displaying main statistics of the SFT fit test

## Examples

``` r
set.seed(12)
data(zma.se)
filt.zma <- filter_by_variance(zma.se, n=500)
zma.set1 <- filt.zma[, sample(colnames(filt.zma), size=20, replace=FALSE)]
zma.set2 <- filt.zma[, sample(colnames(filt.zma), size=20, replace=FALSE)]
list.sets <- list(zma.set1, zma.set2)
cons_sft <- consensus_SFT_fit(list.sets, setLabels = c("Maize1", "Maize2"),
                              cor_method = "pearson")
#>    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
#> 1      5    0.162 -0.241         0.0679    49.8     39.10  120.0
#> 2      6    0.317 -0.353         0.1980    42.9     31.60  109.0
#> 3      7    0.504 -0.441         0.4200    37.7     26.10  101.0
#> 4      8    0.653 -0.529         0.5950    33.5     22.50   93.7
#> 5      9    0.720 -0.598         0.6690    30.2     19.90   87.5
#> 6     10    0.731 -0.656         0.6590    27.4     18.80   82.1
#> 7     11    0.872 -0.686         0.8430    25.1     17.20   77.4
#> 8     12    0.826 -0.721         0.7800    23.1     15.40   73.1
#> 9     13    0.832 -0.739         0.7850    21.4     13.70   69.3
#> 10    14    0.903 -0.737         0.8780    19.9     12.80   65.8
#> 11    15    0.875 -0.756         0.8420    18.6     11.70   62.7
#> 12    16    0.894 -0.764         0.8680    17.4     10.90   59.8
#> 13    17    0.908 -0.756         0.8860    16.4     10.10   57.2
#> 14    18    0.930 -0.761         0.9180    15.4      9.36   54.7
#> 15    19    0.941 -0.762         0.9300    14.6      8.73   52.4
#> 16    20    0.944 -0.764         0.9340    13.8      8.20   50.3
#>    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
#> 1      5   0.0693 -0.178         -0.154    59.7      44.3  137.0
#> 2      6   0.1920 -0.282          0.026    51.9      37.8  127.0
#> 3      7   0.3530 -0.378          0.226    45.8      33.3  119.0
#> 4      8   0.4680 -0.453          0.335    40.9      28.6  111.0
#> 5      9   0.6020 -0.530          0.505    36.9      25.6  105.0
#> 6     10   0.6560 -0.596          0.564    33.6      22.3   98.9
#> 7     11   0.7400 -0.664          0.667    30.8      20.1   93.8
#> 8     12   0.7690 -0.699          0.704    28.3      19.1   89.1
#> 9     13   0.8040 -0.747          0.749    26.2      17.5   84.9
#> 10    14   0.8060 -0.781          0.750    24.4      16.1   81.1
#> 11    15   0.8250 -0.792          0.776    22.7      14.8   77.5
#> 12    16   0.8510 -0.822          0.811    21.3      13.6   74.3
#> 13    17   0.8630 -0.816          0.824    20.0      12.5   71.4
#> 14    18   0.8360 -0.847          0.789    18.8      11.5   68.7
#> 15    19   0.8470 -0.867          0.804    17.7      11.0   66.1
#> 16    20   0.8720 -0.884          0.836    16.8      10.2   63.8
```
