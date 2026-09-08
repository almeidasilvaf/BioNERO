# Filter a gene regulatory network based on optimal scale-free topology fit

Filter a gene regulatory network based on optimal scale-free topology
fit

## Usage

``` r
grn_filter(edgelist, nsplit = 10, bp_param = BiocParallel::SerialParam())
```

## Arguments

- edgelist:

  A gene regulatory network represented as an edge list.

- nsplit:

  Number of groups in which the edge list will be split. Default: 10.

- bp_param:

  BiocParallel back-end to be used. Default: BiocParallel::SerialParam()

## Value

The edge list that best fits the scale-free topology.

## Details

The edge list will be split in n groups and the scale-free topology fit
will be tested for each subset of the edge list. For instance, if an
edge list of 10000 rows is used as input, the function will test SFT fit
for the top 1000 edges, then top 2000 edges, and so on up to the whole
edge list.

## Examples

``` r
data(filt.se)
tfs <- sample(rownames(filt.se), size=50, replace=FALSE)
grn_list <- grn_combined(filt.se, regulators=tfs, nTrees=2)
ranked_grn <- grn_average_rank(grn_list)
# split in only 2 groups for demonstration purposes
filtered_edges <- grn_filter(ranked_grn, nsplit=2)
#> The top number of edges that best fits the scale-free topology is 277
```
