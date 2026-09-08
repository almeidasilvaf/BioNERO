# Get hubs for gene regulatory network

Get hubs for gene regulatory network

Get hubs for protein-protein interaction network

## Usage

``` r
get_hubs_grn(
  edgelist,
  top_percentile = 0.1,
  top_n = NULL,
  return_degree = FALSE,
  ranked = TRUE
)

get_hubs_ppi(
  edgelist,
  top_percentile = 0.1,
  top_n = NULL,
  return_degree = FALSE
)
```

## Arguments

- edgelist:

  A protein-protein interaction network represented as an edge list.

- top_percentile:

  Numeric from 0 to 1 indicating the percentage of proteins with the
  highest degree to consider hubs. Default: 0.1.

- top_n:

  Numeric indicating the number of proteins with the highest degree to
  consider hubs.

- return_degree:

  Logical indicating whether to return a data frame of degree for all
  proteins. If TRUE, the function will return a list instead of a data
  frame. Default: FALSE.

- ranked:

  Logical indicating whether to treat third column of the edge list
  (edge weights) as ranked values. Ignored if the edge list only
  contains 2 columns. Default: TRUE.

## Value

A data frame with gene ID in the first column and out degree in the
second column or a list of two data frames with hubs and degree for all
genes, respectively.

A data frame with protein ID in the first column and degree in the
second column or a list of two data frames with hubs and degree for all
genes, respectively.

## Examples

``` r
data(filt.se)
tfs <- sample(rownames(filt.se), size=50, replace=FALSE)
grn_list <- grn_combined(filt.se, regulators=tfs, nTrees=2)
ranked_grn <- grn_average_rank(grn_list)
# split in only 2 groups for demonstration purposes
filtered_edges <- grn_filter(ranked_grn, nsplit=2)
#> The top number of edges that best fits the scale-free topology is 239
hubs <- get_hubs_grn(filtered_edges)
ppi_edges <- igraph::sample_pa(n = 500)
ppi_edges <- igraph::as_edgelist(ppi_edges)
hubs <- get_hubs_ppi(ppi_edges, return_degree = TRUE)
```
