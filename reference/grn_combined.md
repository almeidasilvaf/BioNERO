# Infer gene regulatory network with multiple algorithms and combine results in a list

Infer gene regulatory network with multiple algorithms and combine
results in a list

## Usage

``` r
grn_combined(
  exp,
  regulators = NULL,
  eps = 0.1,
  estimator_aracne = "spearman",
  estimator_clr = "pearson",
  remove_zero = TRUE,
  ...
)
```

## Arguments

- exp:

  A gene expression data frame with genes in row names and samples in
  column names or a \`SummarizedExperiment\` object.

- regulators:

  A character vector of regulators (e.g., transcription factors or
  miRNAs). All regulators must be included in \`exp\`.

- eps:

  Numeric value indicating the threshold used when removing an edge: for
  each triplet of nodes (i,j,k), the weakest edge, say (ij), is removed
  if its weight is below `min{(ik),(jk)}` - eps. Default: 0.1.

- estimator_aracne:

  Entropy estimator to be used in ARACNE inference. One of
  "mi.empirical", "mi.mm", "mi.shrink", "mi.sg", "pearson", "spearman",
  or "kendall". Default: "spearman".

- estimator_clr:

  Entropy estimator to be used in CLR inference. One of "mi.empirical",
  "mi.mm", "mi.shrink", "mi.sg", "pearson", "spearman", or "kendall".
  Default: "pearson".

- remove_zero:

  Logical indicating whether to remove edges whose weight is exactly
  zero. Zero values indicate edges that were removed by ARACNE. Default:
  TRUE.

- ...:

  Additional arguments passed to \`GENIE3::GENIE3()\`.

## Value

A list of data frames representing edge lists. Each list element is an
edge list for a specific method.

## Examples

``` r
data(filt.se)
tfs <- sample(rownames(filt.se), size=50, replace=FALSE)
grn_list <- grn_combined(filt.se, regulators=tfs, nTrees=2)
```
