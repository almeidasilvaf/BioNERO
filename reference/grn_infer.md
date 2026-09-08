# Infer gene regulatory network with one of three algorithms

The available algorithms are Context Likelihood of Relatedness (CLR),
ARACNE, or GENIE3.

## Usage

``` r
grn_infer(
  exp,
  regulators = NULL,
  method = c("clr", "aracne", "genie3"),
  estimator_clr = "pearson",
  estimator_aracne = "spearman",
  eps = 0.1,
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

- method:

  GRN inference algorithm to be used. One of "clr", "aracne", or
  "genie3".

- estimator_clr:

  Entropy estimator to be used. One of "mi.empirical", "mi.mm",
  "mi.shrink", "mi.sg", "pearson", "spearman", or "kendall". Default:
  "pearson".

- estimator_aracne:

  Entropy estimator to be used. One of "mi.empirical", "mi.mm",
  "mi.shrink", "mi.sg", "pearson", "spearman", or "kendall". Default:
  "spearman".

- eps:

  Numeric value indicating the threshold used when removing an edge: for
  each triplet of nodes (i,j,k), the weakest edge, say (ij), is removed
  if its weight is below `min{(ik),(jk)}` - eps. Default: 0.1.

- remove_zero:

  Logical indicating whether to remove edges whose weight is exactly
  zero. Default: TRUE

- ...:

  Additional arguments passed to \`GENIE3::GENIE3()\`.

## Value

A gene regulatory network represented as an edge list.

## Examples

``` r
data(filt.se)
tfs <- sample(rownames(filt.se), size=20, replace=FALSE)
clr <- grn_infer(filt.se, method = "clr", regulators=tfs)
aracne <- grn_infer(filt.se, method = "aracne", regulators=tfs)
# only 2 trees for demonstration purposes
genie3 <- grn_infer(filt.se, method = "genie3", regulators=tfs, nTrees=2)
```
