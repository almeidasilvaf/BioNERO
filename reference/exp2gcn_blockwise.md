# Infer gene coexpression network from gene expression in a blockwise manner

Infer gene coexpression network from gene expression in a blockwise
manner

## Usage

``` r
exp2gcn_blockwise(
  exp,
  net_type = "signed",
  module_merging_threshold = 0.8,
  SFTpower = NULL,
  cor_method = "pearson",
  TOM_type = NULL,
  max_block_size = 5000,
  min_module_size = 30,
  ...
)
```

## Arguments

- exp:

  Either a \`SummarizedExperiment\` object, or a gene expression
  matrix/data frame with genes in row names and samples in column names.

- net_type:

  Character indicating the type of network to infer. One of 'signed',
  'signed hybrid' or 'unsigned'. Default: 'signed'.

- module_merging_threshold:

  Numeric indicating the minimum correlation threshold to merge similar
  modules into a single one. Default: 0.8.

- SFTpower:

  Numeric scalar indicating the value of the \\\beta\\ power to which
  correlation coefficients will be raised to ensure scale-free topology
  fit. This value can be obtained with the function
  [`SFT_fit()`](SFT_fit.md).

- cor_method:

  Character with correlation method to use. One of "pearson" or
  "biweight". Default: "pearson".

- TOM_type:

  Character specifying the method to use to calculate a topological
  overlap matrix (TOM). If NULL, TOM type will be automatically inferred
  from network type specified in **net_type**. Default: NULL.

- max_block_size:

  Numeric indicating the maximum block size for module detection.

- min_module_size:

  Numeric indicating the minimum module size. Default: 30.

- ...:

  Additional arguments to
  [`WGCNA::blockwiseModules()`](https://rdrr.io/pkg/WGCNA/man/blockwiseModules.html).

## Value

A list containing the following elements:

- *MEs* Data frame of module eigengenes, with samples in rows, and
  module eigengenes in columns.

- *genes_and_modules* Data frame with columns 'Genes' (character) and
  'Modules' (character) indicating the genes and the modules to which
  they belong.

- *params* List with network inference parameters passed as input.

- *dendro_plot_objects* List with objects to plot the dendrogram in
  `plot_dendro_and_colors`. Elements are named 'tree' (an hclust object
  with gene dendrogram), 'Unmerged' (character with per-gene module
  assignments before merging similar modules), and 'Merged' (character
  with per-gene module assignments after merging similar modules).

## Author

Fabricio Almeida-Silva

## Examples

``` r
data(filt.se)
# The SFT fit was previously calculated and the optimal power was 16
cor <- WGCNA::cor
gcn <- exp2gcn_blockwise(
    exp = filt.se, SFTpower = 18, cor_method = "pearson"
)
```
