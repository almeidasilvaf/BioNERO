# Calculate gene significance for a given group of genes

Calculate gene significance for a given group of genes

## Usage

``` r
gene_significance(
  exp,
  metadata,
  metadata_cols = NULL,
  genes = NULL,
  alpha = 0.05,
  cor_method = "pearson",
  min_cor = 0.2,
  use_abs = TRUE
)
```

## Arguments

- exp:

  A gene expression data frame with genes in row names and samples in
  column names or a \`SummarizedExperiment\` object.

- metadata:

  A data frame containing sample names in row names and sample
  annotation in the first column. Ignored if \`exp\` is a
  \`SummarizedExperiment\` object, since the function will extract
  colData.

- metadata_cols:

  A vector (either numeric or character) indicating which columns should
  be extracted from column metadata if **exp** is a
  \`SummarizedExperiment\` object. The vector can contain column indices
  (numeric) or column names (character). By default, all columns are
  used.

- genes:

  Character vector of genes to be correlated with traits. If not given,
  all genes in \`exp\` will be considered.

- alpha:

  Significance level. Default is 0.05.

- cor_method:

  Method to calculate correlation. One of 'pearson', 'spearman' or
  'kendall'. Default is 'spearman'.

- min_cor:

  Minimum correlation coefficient. Default is 0.2.

- use_abs:

  Logical indicating whether to filter by correlation using absolute
  value or not. If TRUE, a `min_cor` of say 0.2 would keep all
  correlations above 0.2 and below -0.2. Default is TRUE.

## Value

A data frame with correlation and correlation p-values for each pair of
gene and trait, with the following variables:

- gene:

  Factor, gene ID.

- trait:

  Factor, trait name. Each trait corresponds to a variable of the sample
  metadata (if numeric) or levels of a variable (if categorical).

- cor:

  Numeric, correlation.

- pvalue:

  Numeric, correlation P-values.

- group:

  Character, name of the metadata variable.

## Author

Fabricio Almeida-Silva

## Examples

``` r
data(filt.se)
gs <- gene_significance(filt.se)
```
