# Collapse gene-level expression data to orthogroup level

For a given list of expression data, this function replaces genes with
their corresponding orthogroups to allow inter-species comparisons.

## Usage

``` r
exp_genes2orthogroups(explist = NULL, og = NULL, summarize = "median")
```

## Arguments

- explist:

  List of expression data frames or SummarizedExperiment objects.

- og:

  Data frame of 3 columns corresponding to orthogroup, species ID, and
  gene ID, respectively. Species IDs must be the same as the names of
  the expression list.

- summarize:

  Centrality measure to summarize multiple paralogous genes in the same
  orthogroup. One of "median" or "mean". Default: "median".

## Value

List of expression data frames for each species with expression
summarized at the orthogroup level.

## Examples

``` r
# \donttest{
data(og.zma.osa)
data(zma.se)
data(osa.se)
explist <- list(zma = zma.se,
                osa = osa.se)
og <- og.zma.osa
exp_ortho <- exp_genes2orthogroups(explist, og, summarize = "mean")
# }
```
