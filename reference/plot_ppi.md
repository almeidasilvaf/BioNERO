# Plot protein-protein interaction network from edge list

Plot protein-protein interaction network from edge list

## Usage

``` r
plot_ppi(
  edgelist_int,
  color_by = "community",
  clustering_method = igraph::cluster_infomap,
  show_labels = "tophubs",
  top_n_hubs = 5,
  add_color_legend = TRUE,
  curvature = 0,
  interactive = FALSE,
  dim_interactive = c(600, 600)
)
```

## Arguments

- edgelist_int:

  Data frame containing the edge list for the PPI network. First column
  is the protein 1 and second column is the protein 2. All other columns
  are interpreted as edge attributes.

- color_by:

  How should nodes be colored? It must be either "community" or a
  2-column data frame containing proteins in the first column and a
  custom annotation in the second column. If "community", a clustering
  algorithm will be applied. Default: "community".

- clustering_method:

  igraph function to be used for community detection. Available
  functions are cluster_infomap, cluster_edge_betweenness,
  cluster_fast_greedy, cluster_walktrap, cluster_spinglass,
  cluster_leading_eigen, cluster_louvain, and cluster_label_prop.
  Default is cluster_infomap.

- show_labels:

  Character indicating which nodes will be labeled. One of "all",
  "allhubs", "tophubs", or "none".

- top_n_hubs:

  Number of top hubs to be labeled. It is only valid if `show_labels`
  equals "tophubs". Default is 5.

- add_color_legend:

  Logical indicating whether to add a color legend for nodes. Default:
  TRUE.

- curvature:

  Numeric indicating the amount of curvature in edges. Negative values
  produce left-hand curves, positive values produce right-hand curves,
  and zero produces a straight line. Default: 0.

- interactive:

  Logical indicating whether the network should be interactive or not.
  Default is FALSE.

- dim_interactive:

  Numeric vector with width and height of window for interactive
  plotting. Default: c(600,600).

## Value

A ggplot object.

## Author

Fabricio Almeida-Silva

## Examples

``` r
ppi_edges <- igraph::sample_pa(n = 500)
ppi_edges <- igraph::as_edgelist(ppi_edges)
p <- plot_ppi(ppi_edges, add_color_legend = FALSE)
```
