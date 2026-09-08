# Plot gene regulatory network from edge list

Plot gene regulatory network from edge list

## Usage

``` r
plot_grn(
  edgelist_grn,
  show_labels = "tophubs",
  top_n_hubs = 5,
  layout = igraph::with_kk,
  arrow.gap = 0.01,
  ranked = TRUE,
  curvature = 0.1,
  interactive = FALSE,
  dim_interactive = c(600, 600)
)
```

## Arguments

- edgelist_grn:

  Data frame containing the edge list for the GRN network. First column
  is the TF and second column is the target gene. All other columns are
  interpreted as edge attributes.

- show_labels:

  Character indicating which nodes will be labeled. One of "all",
  "allhubs", "tophubs", or "none".

- top_n_hubs:

  Number of top hubs to be labeled. It is only valid if `show_labels`
  equals "tophubs". Default is 5.

- layout:

  igraph function for the network layout. One of with_dh, with_drl,
  with_gem, with_lgl, with_fr, with_graphopt, with_kk and with_mds.
  Default is with_kk.

- arrow.gap:

  Numeric indicating the distance between nodes and arrows. Default is
  0.2.

- ranked:

  Logical indicating whether to treat third column of the edge list
  (edge weights) as ranked values. Default: TRUE.

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

A ggplot object containing the network.

## Author

Fabricio Almeida-Silva

## Examples

``` r
data(filt.se)
tfs <- sample(rownames(filt.se), size = 50, replace = FALSE)
grn_edges <- grn_infer(filt.se, method = "clr", regulators = tfs)
p <- plot_grn(grn_edges, ranked = FALSE)
```
