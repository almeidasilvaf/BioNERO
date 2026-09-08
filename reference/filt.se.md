# Filtered maize gene expression data from Shin et al., 2021.

Filtered expression data in transcripts per million (TPM) from Shin et
al., 2021. This is the same data set described in `zma.se`, but it only
contains the top 500 genes with the highest variances. This data set was
created to be used in unit tests and examples.

## Usage

``` r
data(filt.se)
```

## Format

An object of class `SummarizedExperiment`

## References

Shin, J., Marx, H., Richards, A., Vaneechoutte, D., Jayaraman, D.,
Maeda, J., ... & Roy, S. (2021). A network-based comparative framework
to study conservation and divergence of proteomes in plant phylogenies.
Nucleic Acids Research, 49(1), e3-e3.

## Examples

``` r
data(filt.se)
```
