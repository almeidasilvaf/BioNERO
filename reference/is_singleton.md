# Logical expression to check if gene or gene set is singleton or not

Logical expression to check if gene or gene set is singleton or not

## Usage

``` r
is_singleton(genes, og)
```

## Arguments

- genes:

  Character containing gene or group of genes to be evaluated.

- og:

  Data frame of 3 columns corresponding to orthogroup, species ID, and
  gene ID, respectively.

## Value

Vector of logical values indicating if gene or group of genes is
singleton or not.

## See also

`is_duplicated`

## Author

Fabricio Almeida-Silva

## Examples

``` r
data(og.zma.osa)
data(filt.se)
genes <- tail(rownames(filt.se), n = 100)
is_singleton(genes, og.zma.osa)
#> Zm00001d042050 Zm00001d042084 Zm00001d042308 Zm00001d042453 Zm00001d042525 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d042730 Zm00001d042772 Zm00001d042935 Zm00001d042966 Zm00001d043170 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d043299 Zm00001d043382 Zm00001d043465 Zm00001d043606 Zm00001d043942 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d043998 Zm00001d044130 Zm00001d044228 Zm00001d044246 Zm00001d044287 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d044664 Zm00001d044685 Zm00001d044686 Zm00001d044747 Zm00001d045000 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d045025 Zm00001d045042 Zm00001d045139 Zm00001d045323 Zm00001d045431 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d045448 Zm00001d045544 Zm00001d045888 Zm00001d046378 Zm00001d046449 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d046555 Zm00001d046583 Zm00001d046767 Zm00001d047203 Zm00001d047253 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d047282 Zm00001d047296 Zm00001d047479 Zm00001d047581 Zm00001d047697 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d047707 Zm00001d047765 Zm00001d047787 Zm00001d048091 Zm00001d048201 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d048346 Zm00001d048368 Zm00001d048491 Zm00001d048611 Zm00001d048787 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d049099 Zm00001d049166 Zm00001d049239 Zm00001d049500 Zm00001d049541 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d049552 Zm00001d049641 Zm00001d049674 Zm00001d049790 Zm00001d049826 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d049933 Zm00001d050032 Zm00001d050100 Zm00001d050193 Zm00001d050375 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d050961 Zm00001d051001 Zm00001d051056 Zm00001d051420 Zm00001d051478 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d051543 Zm00001d051591 Zm00001d051595 Zm00001d051659 Zm00001d051685 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d052101 Zm00001d052216 Zm00001d052718 Zm00001d052720 Zm00001d052855 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d052885 Zm00001d053090 Zm00001d053103 Zm00001d053228 Zm00001d053239 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d053354 Zm00001d053593 Zm00001d053632 Zm00001d053633 Zm00001d053635 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> Zm00001d053636 Zm00001d053834 Zm00001d053838 Zm00001d053845 Zm00001d053863 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
```
