# Preprocess expression data for network reconstruction

Preprocess expression data for network reconstruction

## Usage

``` r
exp_preprocess(
  exp,
  NA_rm = TRUE,
  replaceby = 0,
  Zk_filtering = TRUE,
  zk = -2,
  cor_method = "spearman",
  remove_nonexpressed = TRUE,
  method = "median",
  min_exp = 1,
  min_percentage_samples = 0.25,
  remove_confounders = TRUE,
  variance_filter = FALSE,
  n = NULL,
  percentile = NULL,
  vstransform = FALSE
)
```

## Arguments

- exp:

  A gene expression data frame with genes in row names and samples in
  column names or a \`SummarizedExperiment\` object.

- NA_rm:

  Logical. It specifies whether to remove missing values from the
  expression data frame or not. Default = TRUE.

- replaceby:

  If NA_rm is TRUE, what to use instead of NAs. One of 0 or 'mean'.
  Default is 0.

- Zk_filtering:

  Logical. It specifies whether to filter outlying samples by Zk or not.
  Default: TRUE.

- zk:

  If Zk_filtering is TRUE, the standardized connectivity threshold.
  Samples below this threshold will be considered outliers. Default is
  -2.

- cor_method:

  If Zk_filtering is TRUE, the correlation method to use. One of
  'spearman', 'bicor', or 'pearson'. Default is 'spearman'.

- remove_nonexpressed:

  Logical. It specifies whether non-expressed genes should be removed or
  not. Default is TRUE.

- method:

  If remove_nonexpressed is TRUE, the criterion to filter non-expressed
  genes out. One of "mean", "median", "percentage", or "allsamples".
  Default is 'median'.

- min_exp:

  If method is 'mean', 'median', or 'allsamples', the minimum value for
  a gene to be considered expressed. If method is 'percentage', the
  minimum value each gene must have in at least n percent of samples to
  be considered expressed.

- min_percentage_samples:

  If method is 'percentage', expressed genes must have expression \>=
  min_exp in at least this percentage. Values must range from 0 to 1.
  Default = 0.25.

- remove_confounders:

  Logical. If TRUE, it removes principal components that add noise to
  the data.

- variance_filter:

  Logical. If TRUE, it will filter genes by variance. Default is FALSE.

- n:

  If variance_filter is TRUE, the number of most variable genes to keep.

- percentile:

  If variance_filter is TRUE, the percentage of most variable genes to
  keep.

- vstransform:

  Logical indicating if data should be variance stabilizing transformed.
  This parameter can only be set to TRUE if data is a matrix of raw read
  counts.

## Value

Processed gene expression data frame with gene IDs in row names and
sample names in column names or \`SummarizedExperiment\` object.

## References

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of
fold change and dispersion for RNA-seq data with DESeq2. Genome biology,
15(12), 1-21.

## See also

[`varianceStabilizingTransformation`](https://rdrr.io/pkg/DESeq2/man/varianceStabilizingTransformation.html)

## Author

Fabricio Almeida-Silva

## Examples

``` r
data(zma.se)
exp <- exp_preprocess(zma.se, variance_filter=TRUE, n=1000)
#> Number of removed samples: 1
```
