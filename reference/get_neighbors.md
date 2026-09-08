# Get 1st-order neighbors of a given gene or group of genes

Get 1st-order neighbors of a given gene or group of genes

## Usage

``` r
get_neighbors(genes, net, cor_threshold = 0.7)
```

## Arguments

- genes:

  Character vector containing genes from which direct neighbors will be
  extracted.

- net:

  List object returned by `exp2gcn`.

- cor_threshold:

  Correlation threshold to filter connections. As a weighted network is
  a fully connected graph, a cutoff must be selected. Default is 0.7.

## Value

List containing 1st-order neighbors for each input gene.

## See also

`exp2gcn` `SFT_fit`

## Author

Fabricio Almeida-Silva

## Examples

``` r
data(filt.se)
genes <- rownames(filt.se)[1:10]
gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson")
#> ..connectivity..
#> ..matrix multiplication (system BLAS)..
#> ..normalization..
#> ..done.
neighbors <- get_neighbors(genes, gcn)
```
