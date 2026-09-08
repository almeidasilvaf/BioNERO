# Detect communities in a network

Detect communities in a network

## Usage

``` r
detect_communities(edgelist, method = igraph::cluster_infomap, directed = TRUE)
```

## Arguments

- edgelist:

  Data frame containing the network as an edge list. First column must
  be node 1 and second column must be node 2. Additional columns will be
  interpreted as edge attributes and will be modified by this function.

- method:

  igraph function to be used for community detection. Available
  functions are cluster_infomap, cluster_edge_betweenness,
  cluster_fast_greedy, cluster_walktrap, cluster_spinglass,
  cluster_leading_eigen, cluster_louvain, and cluster_label_prop.
  Default is cluster_infomap.

- directed:

  Logical indicating whether the network is directed (GRN only) or not
  (GCN and PPI networks). Default: TRUE.

## Value

A data frame containing node names in the first column, and communities
to which nodes belong in the second column.

## See also

[`cluster_infomap`](https://r.igraph.org/reference/cluster_infomap.html),
[`cluster_edge_betweenness`](https://r.igraph.org/reference/cluster_edge_betweenness.html),
[`cluster_fast_greedy`](https://r.igraph.org/reference/cluster_fast_greedy.html),
[`cluster_walktrap`](https://r.igraph.org/reference/cluster_walktrap.html),
[`cluster_spinglass`](https://r.igraph.org/reference/cluster_spinglass.html),
[`cluster_leading_eigen`](https://r.igraph.org/reference/cluster_leading_eigen.html),
[`cluster_louvain`](https://r.igraph.org/reference/cluster_louvain.html),
[`cluster_label_prop`](https://r.igraph.org/reference/cluster_label_prop.html)

## Author

Fabricio Almeida-Silva

## Examples

``` r
data(filt.se)
tfs <- sample(rownames(filt.se), size=50, replace=FALSE)
grn_edges <- grn_infer(filt.se, method = "clr", regulators = tfs)
com <- detect_communities(grn_edges, directed=TRUE)
```
