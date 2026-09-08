# Calculate network statistics

Calculate network statistics

## Usage

``` r
net_stats(
  adj_matrix = NULL,
  net_type = c("gcn", "ppi", "grn"),
  calculate_additional = FALSE
)
```

## Arguments

- adj_matrix:

  Adjacency matrix that represents the network.

- net_type:

  One of "gcn" (gene coexpression network), "ppi" (protein-protein
  interaction), or "grn" (gene regulatory network).

- calculate_additional:

  Logical indicating whether to calculate additional network statistics
  (betweenness and closeness). Default is FALSE.

## Value

A list containing the following elements:

- Connectivity

- ScaledConnectivity

- ClusterCoef

- MAR (for gcn only)

- Density

- Centralization

- Heterogeneity (gcn only)

- Diameter

- Betweenness

- Closeness

## See also

[`graph_from_adjacency_matrix`](https://r.igraph.org/reference/graph_from_adjacency_matrix.html),
[`cliques`](https://r.igraph.org/reference/cliques.html),[`diameter`](https://r.igraph.org/reference/diameter.html),
[`estimate_betweenness`](https://r.igraph.org/reference/estimate_betweenness.html),[`V`](https://r.igraph.org/reference/V.html),
[`closeness`](https://r.igraph.org/reference/closeness.html),[`degree`](https://r.igraph.org/reference/degree.html),
[`transitivity`](https://r.igraph.org/reference/transitivity.html),[`edge_density`](https://r.igraph.org/reference/edge_density.html),
[`centr_degree`](https://r.igraph.org/reference/centr_degree.html)
[`fundamentalNetworkConcepts`](https://rdrr.io/pkg/WGCNA/man/fundamentalNetworkConcepts.html)

## Examples

``` r
# \donttest{
data(filt.se)
set.seed(12)
filt.se <- exp_preprocess(
    filt.se, Zk_filtering = FALSE, variance_filter = TRUE, n = 200
)
gcn <- exp2gcn(
    filt.se, SFTpower = 7, cor_method = "pearson", net_type = "signed hybrid"
)
#> ..connectivity..
#> ..matrix multiplication (system BLAS)..
#> ..normalization..
#> ..done.
stats <- net_stats(gcn$adjacency_matrix, net_type = "gcn")
# }
```
