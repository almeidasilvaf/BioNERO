# Identify consensus modules across independent data sets

Identify consensus modules across independent data sets

## Usage

``` r
consensus_modules(
  exp_list,
  metadata,
  power,
  cor_method = "spearman",
  net_type = "signed hybrid",
  module_merging_threshold = 0.8,
  TOM_type = NULL,
  verbose = FALSE
)
```

## Arguments

- exp_list:

  A list containing the expression data frames with genes in row names
  and samples in column names or \`SummarizedExperiment\` objects. The
  list can be created by using `list(exp1, exp2, ..., expn)`.

- metadata:

  A data frame containing sample names in row names and sample
  annotation in the first column. Ignored if \`exp_list\` is a list of
  \`SummarizedExperiment\` objects, since the function will extract
  colData.

- power:

  Numeric vector of beta power for each expression set as calculated by
  `consensus_SFT_fit`.

- cor_method:

  Correlation method used for network reconstruction. One of "spearman"
  (default), "biweight", or "pearson".

- net_type:

  Network type. One of "signed hybrid" (default), "signed" or
  "unsigned".

- module_merging_threshold:

  Correlation threshold to merge similar modules into a single one.
  Default: 0.8.

- TOM_type:

  Character indicating the type of Topological Overlap Matrix to (TOM)
  create. One of 'unsigned', 'signed', 'signed Nowick', 'unsigned 2',
  'signed 2', and 'signed Nowick 2'. By default, TOM type is
  automatically selected based on network type.

- verbose:

  Logical indicating whether to display progress messages or not.
  Default: FALSE.

## Value

A list containing 4 elements:

- consMEs:

  Consensus module eigengenes

- exprSize:

  Description of the multi-set object returned by the function
  [`WGCNA::checkSets`](https://rdrr.io/pkg/WGCNA/man/checkSets.html)

- sampleInfo:

  Metadata for each expression set

- genes_cmodules:

  Data frame of genes and consensus modules

- dendro_plot_objects:

  Objects to be used in dendrogram plotting

## Examples

``` r
set.seed(12)
data(zma.se)
filt.zma <- filter_by_variance(zma.se, n=500)
zma.set1 <- filt.zma[, sample(colnames(filt.zma), size=20, replace=FALSE)]
zma.set2 <- filt.zma[, sample(colnames(filt.zma), size=20, replace=FALSE)]
list.sets <- list(zma.set1, zma.set2)
# SFT power previously identified with consensus_SFT_fit()
cons_mod <- consensus_modules(list.sets, power = c(11, 13),
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
```
