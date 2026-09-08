# Rank edge weights for GRNs and calculate average across different methods

Rank edge weights for GRNs and calculate average across different
methods

## Usage

``` r
grn_average_rank(list_edges)
```

## Arguments

- list_edges:

  List containing edge lists as returned by the function `grn_combined`.

## Value

Edge list containing regulator, target and mean rank from all
algorithms.

## Examples

``` r
data(filt.se)
tfs <- sample(rownames(filt.se), size=50, replace=FALSE)
grn_list <- grn_combined(filt.se, regulators=tfs, nTrees=2)
ranked_grn <- grn_average_rank(grn_list)
```
