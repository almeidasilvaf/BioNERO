# Get edge list from an adjacency matrix for a group of genes

Get edge list from an adjacency matrix for a group of genes

## Usage

``` r
get_edge_list(
  net,
  genes = NULL,
  module = NULL,
  filter = FALSE,
  method = "optimalSFT",
  r_optimal_test = seq(0.4, 0.9, by = 0.1),
  Zcutoff = 1.96,
  pvalue_cutoff = 0.05,
  rcutoff = 0.7,
  nSamples = NULL,
  check_SFT = FALSE,
  bp_param = BiocParallel::SerialParam()
)
```

## Arguments

- net:

  List object returned by `exp2gcn`.

- genes:

  Character vector containing a subset of genes from which edges will be
  extracted. It can be ignored if the user wants to extract an edge list
  for a given module instead of individual genes.

- module:

  Character with module name from which edges will be extracted. To
  include 2 or more modules, input the names in a character vector.

- filter:

  Logical indicating whether to filter the edge list or not.

- method:

  Method to filter spurious correlations. One of "Zscore", "optimalSFT",
  "pvalue" or "min_cor". See details for more information on the
  methods. Default: 'optimalSFT'

- r_optimal_test:

  Numeric vector with the correlation thresholds to be tested for
  optimal scale-free topology fit. Only valid if `method` equals
  "optimalSFT". Default: seq(0.4, 0.9, by = 0.1)

- Zcutoff:

  Minimum Z-score threshold. Only valid if `method` equals "Zscore".
  Default: 1.96

- pvalue_cutoff:

  Maximum P-value threshold. Only valid if `method` equals "pvalue".
  Default: 0.05

- rcutoff:

  Minimum correlation threshold. Only valid if `method` equals
  "min_cor". Default: 0.7

- nSamples:

  Number of samples in the data set from which the correlation matrix
  was calculated. Only required if `method` equals "pvalue".

- check_SFT:

  Logical indicating whether to test if the resulting network is close
  to a scale-free topology or not. Default: FALSE.

- bp_param:

  BiocParallel back-end to be used. Default: BiocParallel::SerialParam()

## Value

Data frame with edge list for the input genes.

## Details

The default method ("optimalSFT") will create several different edge
lists by filtering the original correlation matrix by the thresholds
specified in `r_optimal_test`. Then, it will calculate a scale-free
topology fit index for each of the possible networks and return the
network that best fits the scale-free topology. The method "Zscore" will
apply a Fisher Z-transformation for the correlation coefficients and
remove the Z-scores below the threshold specified in `Zcutoff`. The
method "pvalue" will calculate Student asymptotic p-value for the
correlations and remove correlations whose p-values are above the
threshold specified in `pvalue_cutoff`. The method "min_cor" will remove
correlations below the minimum correlation threshold specified in
`rcutoff`.

## See also

[`scaleFreeFitIndex`](https://rdrr.io/pkg/WGCNA/man/scaleFreeFitIndex.html)

`SFT_fit`

`exp2gcn`.

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
genes <- rownames(filt.se)[1:50]
edges <- get_edge_list(gcn, genes=genes, filter = FALSE)
```
