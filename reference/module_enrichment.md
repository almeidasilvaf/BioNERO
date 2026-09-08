# Perform enrichment analysis for coexpression network modules

Perform enrichment analysis for coexpression network modules

## Usage

``` r
module_enrichment(
  net = NULL,
  background_genes,
  annotation,
  column = NULL,
  correction = "BH",
  p = 0.05,
  min_setsize = 10,
  max_setsize = 500,
  bp_param = BiocParallel::SerialParam()
)
```

## Arguments

- net:

  List object returned by `exp2gcn`.

- background_genes:

  Character vector of genes to be used as background for the Fisher's
  Exact Test.

- annotation:

  Annotation data frame with genes in the first column and functional
  annotation in the other columns. This data frame can be exported from
  Biomart or similar databases.

- column:

  Column or columns of `annotation` to be used for enrichment. Both
  character or numeric values with column indices can be used. If users
  want to supply more than one column, input a character or numeric
  vector. Default: all columns from `annotation`.

- correction:

  Multiple testing correction method. One of "holm", "hochberg",
  "hommel", "bonferroni", "BH", "BY", "fdr" or "none". Default is "BH".

- p:

  P-value threshold. P-values below this threshold will be considered
  significant. Default is 0.05.

- min_setsize:

  Numeric indicating the minimum gene set size to be considered. Gene
  sets correspond to levels of each variable in **annotation**).
  Default: 10.

- max_setsize:

  Numeric indicating the maximum gene set size to be considered. Gene
  sets correspond to levels of each variable in **annotation**).
  Default: 500.

- bp_param:

  BiocParallel back-end to be used. Default: BiocParallel::SerialParam()

## Value

A data frame of overrepresentation results with the following variables:

- term:

  character, functional term ID/name.

- genes:

  numeric, intersection length between input genes and genes in a
  particular functional term.

- all:

  numeric, number of all genes in a particular functional term.

- pval:

  numeric, P-value for the hypergeometric test.

- padj:

  numeric, P-value adjusted for multiple comparisons using the method
  specified in parameter **adj**.

- category:

  character, name of the grouping variable (i.e., column name of
  **annotation**).

- module:

  character, module name.

## Author

Fabricio Almeida-Silva

## Examples

``` r
# \donttest{
data(filt.se)
data(zma.interpro)
background <- rownames(filt.se)
gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson")
#> ..connectivity..
#> ..matrix multiplication (system BLAS)..
#> ..normalization..
#> ..done.
mod_enrich <- module_enrichment(gcn, background, zma.interpro, p=1)
#> Enrichment analysis for module black...
#> Enrichment analysis for module blue...
#> Enrichment analysis for module brown...
#> Enrichment analysis for module green...
#> Enrichment analysis for module red...
#> Enrichment analysis for module yellow...
# }
```
