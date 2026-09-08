# Apply Principal Component (PC)-based correction for confounding artifacts

Apply Principal Component (PC)-based correction for confounding
artifacts

## Usage

``` r
PC_correction(exp, verbose = FALSE)
```

## Arguments

- exp:

  A gene expression data frame with genes in row names and samples in
  column names or a \`SummarizedExperiment\` object.

- verbose:

  Logical indicating whether to display progress messages or not.
  Default: FALSE.

## Value

Corrected expression data frame or \`SummarizedExperiment\` object.

## References

Parsana, P., Ruberman, C., Jaffe, A. E., Schatz, M. C., Battle, A., &
Leek, J. T. (2019). Addressing confounding artifacts in reconstruction
of gene co-expression networks. Genome biology, 20(1), 1-6.

## See also

[`num.sv`](https://rdrr.io/pkg/sva/man/num.sv.html),[`sva_network`](https://rdrr.io/pkg/sva/man/sva_network.html)

## Author

Fabricio Almeida-Silva

## Examples

``` r
data(zma.se)
exp <- filter_by_variance(zma.se, n=500)
exp <- PC_correction(exp)
```
