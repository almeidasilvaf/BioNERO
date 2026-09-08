# Calculate pairwise correlations between genes in a matrix

Calculate pairwise correlations between genes in a matrix

## Usage

``` r
exp2cor(exp, cor_method = "pearson")
```

## Arguments

- exp:

  A numeric matrix containing a gene expression matrix, with genes in
  rows and samples in columns.

- cor_method:

  Character indicating the correlation method to use. One of "pearson",
  "spearman", or "biweight". Default: "pearson".

## Value

A numeric, symmetric matrix with pairwise correlations between genes.

## Examples

``` r
# Simulate an expression matrix with 100 genes and 50 samples
exp <- matrix(
    rnorm(100 * 50, mean = 10, sd = 2),
    nrow = 100,
    dimnames = list(
        paste0("gene", seq_len(100)),
        paste0("sample", seq_len(50))
    )
)

# Calculate correlation matrix
cor_mat <- exp2cor(exp)
```
