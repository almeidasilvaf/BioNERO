# Correlate set-specific modules and consensus modules to sample information

Correlate set-specific modules and consensus modules to sample
information

## Usage

``` r
consensus_trait_cor(consensus, cor_method = "pearson", metadata_cols = NULL)
```

## Arguments

- consensus:

  Consensus network returned by `consensus_modules`.

- cor_method:

  Correlation method to be used. One of 'spearman' or 'pearson'.
  Default: 'pearson'.

- metadata_cols:

  A vector (either numeric or character) indicating which columns should
  be extracted from column metadata if **exp** is a
  \`SummarizedExperiment\` object. The vector can contain column indices
  (numeric) or column names (character). By default, all columns are
  used.

## Value

Data frame of consensus module-trait correlations and p-values, with the
following variables:

- trait:

  Factor, trait name. Each trait corresponds to a variable of the sample
  metadata (if numeric) or levels of a variable (if categorical).

- ME:

  Factor, module eigengene.

- cor:

  Numeric, correlation.

- pvalue:

  Numeric, correlation P-values.

- group:

  Character, name of the metadata variable.

## Examples

``` r
set.seed(12)
data(zma.se)
filt.zma <- filter_by_variance(zma.se, n=500)
zma.set1 <- filt.zma[, sample(colnames(filt.zma), size=20, replace=FALSE)]
zma.set2 <- filt.zma[, sample(colnames(filt.zma), size=20, replace=FALSE)]
list.sets <- list(zma.set1, zma.set2)
# SFT power previously identified with consensus_SFT_fit()
consensus <- consensus_modules(list.sets, power = c(11, 13),
                               cor_method = "pearson")
#> ..connectivity..
#> ..matrix multiplication (system BLAS)..
#> ..normalization..
#> ..done.
#> ..connectivity..
#> ..matrix multiplication (system BLAS)..
#> ..normalization..
#> ..done.
#>  ..done.
#>  multiSetMEs: Calculating module MEs.
#>    Working on set 1 ...
#>    Working on set 2 ...
consensus_trait <- consensus_trait_cor(consensus, cor_method = "pearson")
```
