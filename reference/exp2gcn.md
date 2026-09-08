# Infer gene coexpression network from gene expression

Infer gene coexpression network from gene expression

## Usage

``` r
exp2gcn(
  exp,
  net_type = "signed",
  module_merging_threshold = 0.8,
  SFTpower = NULL,
  cor_method = "spearman",
  TOM_type = NULL,
  min_module_size = 30,
  return_cormat = TRUE,
  verbose = FALSE
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

  Character with correlation method to use. One of "pearson", "biweight"
  or "spearman". Default: "spearman".

- TOM_type:

  Character specifying the method to use to calculate a topological
  overlap matrix (TOM). If NULL, TOM type will be automatically inferred
  from network type specified in **net_type**. Default: NULL.

- min_module_size:

  Numeric indicating the minimum module size. Default: 30.

- return_cormat:

  Logical indicating whether the correlation matrix should be returned.
  If TRUE (default), an element named \`correlation_matrix\` containing
  the correlation matrix will be included in the result list.

- verbose:

  Logical indicating whether to display progress messages or not.
  Default: FALSE.

## Value

A list containing the following elements:

- *adjacency_matrix* Numeric matrix with network adjacencies.

- *MEs* Data frame of module eigengenes, with samples in rows, and
  module eigengenes in columns.

- *genes_and_modules* Data frame with columns 'Genes' (character) and
  'Modules' (character) indicating the genes and the modules to which
  they belong.

- *kIN* Data frame of degree centrality for each gene, with columns
  'kTotal' (total degree), 'kWithin' (intramodular degree), 'kOut'
  (extra-modular degree), and 'kDiff' (difference between the intra- and
  extra-modular degree).

- *correlation_matrix* Numeric matrix with pairwise correlation
  coefficients between genes. If parameter **return_cormat** is FALSE,
  this will be NULL.

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
gcn <- exp2gcn(exp = filt.se, SFTpower = 18, cor_method = "pearson")
#> ..connectivity..
#> ..matrix multiplication (system BLAS)..
#> ..normalization..
#> ..done.
```
