# Quantile normalize the expression data

Quantile normalize the expression data

## Usage

``` r
q_normalize(exp)
```

## Arguments

- exp:

  A gene expression data frame with genes in row names and samples in
  column names.

## Value

Expression matrix with normalized values

## Examples

``` r
data(zma.se)
exp <- SummarizedExperiment::assay(zma.se)
norm_exp <- q_normalize(exp)
```
