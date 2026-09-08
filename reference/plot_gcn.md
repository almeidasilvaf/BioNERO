# Plot gene coexpression network from edge list

Plot gene coexpression network from edge list

## Usage

``` r
plot_gcn(
  edgelist_gcn,
  net,
  color_by = "module",
  hubs = NULL,
  show_labels = "tophubs",
  top_n_hubs = 5,
  curvature = 0,
  interactive = FALSE,
  dim_interactive = c(600, 600)
)
```

## Arguments

- edgelist_gcn:

  Data frame containing the edge list for the GCN. The edge list can be
  generated with [`get_edge_list()`](get_edge_list.md).

- net:

  List object returned by `exp2net`.

- color_by:

  How should nodes be colored? It must be either "module" (nodes will
  have the colors of their modules) or a 2-column data frame containing
  genes in the first column and a custom gene annotation in the second
  column. Default: "module".

- hubs:

  Data frame containing hub genes in the first column, their modules in
  the second column, and intramodular connectivity in the third column.

- show_labels:

  Character indicating which nodes will be labeled. One of "all",
  "allhubs", "tophubs", or "none". Default: tophubs.

- top_n_hubs:

  Number of top hubs to be labeled. It is only valid if `show_labels`
  equals "tophubs". Default is 5.

- curvature:

  Numeric indicating the amount of curvature in edges. Negative values
  produce left-hand curves, positive values produce right-hand curves,
  and zero produces a straight line. Default: 0.1.

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
data(filt.se)
gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson")
#> ..connectivity..
#> ..matrix multiplication (system BLAS)..
#> ..normalization..
#> ..done.
gcn_edges <- get_edge_list(gcn, module="brown", filter=TRUE,
                           method="min_cor")
hubs <- get_hubs_gcn(filt.se, gcn)
p <- plot_gcn(gcn_edges, gcn, hubs = hubs)
```
