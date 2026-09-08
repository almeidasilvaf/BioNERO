# Infer gene regulatory network from expression data

Infer gene regulatory network from expression data

## Usage

``` r
exp2grn(
  exp,
  regulators = NULL,
  eps = 0,
  estimator_aracne = "spearman",
  estimator_clr = "pearson",
  remove_zero = TRUE,
  nsplit = 10,
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
  if its weight is below `min{(ik),(jk)}` - eps. Default: 0.

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

- nsplit:

  Number of groups in which the edge list will be split. Default: 10.

- ...:

  Additional arguments passed to \`GENIE3::GENIE3()\`.

## Value

A filtered edge list with regulators in the first column and targets in
the second column.

## Details

This function infers GRNs with ARACNE, GENIE3 and CLR, ranks correlation
weights for each GRN and calculates the average rank for each edge.
Then, the resulting GRN is filtered to keep the top n edges that lead to
the optimal scale-free topology fit.

## Examples

``` r
data(filt.se)
tfs <- sample(rownames(filt.se), size=50, replace=FALSE)
# Test with small number of trees for demonstration purpose
grn <- exp2grn(filt.se, regulators = tfs, nTrees=2, nsplit=2)
#> The top number of edges that best fits the scale-free topology is 254
```
