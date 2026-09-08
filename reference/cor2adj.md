# Calculate an adjacency matrix from a correlation matrix

Calculate an adjacency matrix from a correlation matrix

## Usage

``` r
cor2adj(cor_matrix, beta, net_type = "signed hybrid")
```

## Arguments

- cor_matrix:

  A numeric, symmetric matrix with pairwise correlations between genes
  (i.e., a 'correlation matrix').

- beta:

  Numeric scalar indicating the value of the \\\beta\\ power to which
  correlation coefficients will be raised to ensure scale-free topology
  fit.

- net_type:

  Character indicating the type of network to infer. Default: "signed
  hybrid".

## Value

A numeric, symmetric matrix with network adjacency values between genes.

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

# Calculate adjacency matrix (random value for beta)
adj <- cor2adj(cor_mat, beta = 9)
```
