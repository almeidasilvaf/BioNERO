# Get housekeeping genes from global expression profile

Get housekeeping genes from global expression profile

## Usage

``` r
get_HK(exp)
```

## Arguments

- exp:

  A gene expression data frame with genes in row names and samples in
  column names or a \`SummarizedExperiment\` object.

## Value

Character vector of housekeeping gene IDs.

## Details

This function identifies housekeeping genes, which are broadly expressed
genes with low variation in a global scale across samples. For some
cases, users would want to remove these genes as they are not
interesting for coexpression network analyses. See references for more
details.

## References

Machado, F.B., Moharana, K.C., Almeida‐Silva, F., Gazara, R.K.,
Pedrosa‐Silva, F., Coelho, F.S., Grativol, C. and Venancio, T.M. (2020),
Systematic analysis of 1298 RNA‐Seq samples and construction of a
comprehensive soybean (Glycine max) expression atlas. Plant J, 103:
1894-1909.

## Author

Fabricio Almeida-Silva

## Examples

``` r
data(zma.se)
hk <- get_HK(zma.se)
```
