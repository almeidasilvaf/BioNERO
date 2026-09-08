# Transform a correlation matrix to an edge list

Transform a correlation matrix to an edge list

## Usage

``` r
cormat_to_edgelist(matrix)
```

## Arguments

- matrix:

  Symmetrical correlation matrix.

## Value

A 2-column data frame containing node 1, node 2 and edge weight.

## Examples

``` r
data(filt.se)
cor_mat <- cor(t(SummarizedExperiment::assay(filt.se)))
edgelist <- cormat_to_edgelist(cor_mat)
```
