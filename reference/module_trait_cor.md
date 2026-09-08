# Correlate module eigengenes to trait

Correlate module eigengenes to trait

## Usage

``` r
module_trait_cor(
  exp,
  metadata,
  MEs,
  metadata_cols = NULL,
  cor_method = "pearson"
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

- MEs:

  Module eigengenes. It is the 2nd element of the result list generated
  by the function `exp2gcn`.

- metadata_cols:

  A vector (either numeric or character) indicating which columns should
  be extracted from column metadata if **exp** is a
  \`SummarizedExperiment\` object. The vector can contain column indices
  (numeric) or column names (character). By default, all columns are
  used.

- cor_method:

  Method to calculate correlation. One of 'pearson', 'spearman' or
  'kendall'. Default is 'spearman'.

## Value

A data frame with correlation and correlation p-values for each pair of
ME and trait, with the following variables:

- ME:

  Factor, module eigengene.

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
gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson")
#> ..connectivity..
#> ..matrix multiplication (system BLAS)..
#> ..normalization..
#> ..done.
module_trait_cor(filt.se, MEs = gcn$MEs)
#>          ME          trait         cor     pvalue  group
#> 1   MEblack      endosperm  0.02815842 0.88689440 Tissue
#> 2   MEblack         pollen  0.22278388 0.25449144 Tissue
#> 3   MEblack whole_seedling -0.19659286 0.31601910 Tissue
#> 4    MEblue      endosperm  0.12177943 0.53702286 Tissue
#> 5    MEblue         pollen  0.25484034 0.19062230 Tissue
#> 6    MEblue whole_seedling -0.30601447 0.11325900 Tissue
#> 7   MEbrown      endosperm  0.43581294 0.02043727 Tissue
#> 8   MEbrown         pollen -0.10954013 0.57897711 Tissue
#> 9   MEbrown whole_seedling -0.31064770 0.10762859 Tissue
#> 10  MEgreen      endosperm -0.26401402 0.17460060 Tissue
#> 11  MEgreen         pollen -0.19109098 0.33002030 Tissue
#> 12  MEgreen whole_seedling  0.38589748 0.04253894 Tissue
#> 13   MEgrey      endosperm -0.25127399 0.19711727 Tissue
#> 14   MEgrey         pollen -0.12106772 0.53942247 Tissue
#> 15   MEgrey whole_seedling  0.32058309 0.09626135 Tissue
#> 16    MEred      endosperm -0.09629821 0.62592835 Tissue
#> 17    MEred         pollen  0.19822084 0.31194768 Tissue
#> 18    MEred whole_seedling -0.06499422 0.74246952 Tissue
#> 19 MEyellow      endosperm  0.16325655 0.40649519 Tissue
#> 20 MEyellow         pollen  0.16243302 0.40889882 Tissue
#> 21 MEyellow whole_seedling -0.27262137 0.16045070 Tissue
```
