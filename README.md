
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BioNERO <img src='man/figures/logo.png' align="right" height="139" />

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/almeidasilvaf/BioNERO)](https://github.com/almeidasilvaf/BioNERO/issues)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check-bioc](https://github.com/almeidasilvaf/BioNERO/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/almeidasilvaf/BioNERO/actions)
[![Codecov test
coverage](https://codecov.io/gh/almeidasilvaf/BioNERO/branch/master/graph/badge.svg)](https://codecov.io/gh/almeidasilvaf/BioNERO?branch=master)
<!-- badges: end -->

`BioNERO` aims to integrate all aspects of biological network inference
in a single package, so users don’t have to learn the syntaxes of
several packages and how to communicate among them. `BioNERO` features:

-   **Expression data preprocessing** using state-of-the-art techniques
    for network inference.
-   **Automated exploratory data analyses**, including principal
    component analysis (PCA) and heatmaps of gene expression or sample
    correlations.
-   **Inference of gene coexpression networks (GCNs)** using the popular
    WGCNA algorithm.
-   **Inference of gene regulatory networks (GRNs)** based on the
    “wisdom of the crowds” principle. This principle consists in
    inferring GRNs with multiple algorithms (here, CLR, GENIE3 and
    ARACNE) and calculating the average rank for each interaction pair.
-   **Exploration of network topology** of GCNs, GRNs, and
    protein-protein interaction networks.
-   **Network visualization**.
-   **Network comparison**, including identification of consensus
    modules across independent expression sets, and calculation of intra
    and interspecies module preservation statistics between different
    networks.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `BioNERO` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("BioNERO")
```

And the development version from
[GitHub](https://github.com/almeidasilvaf/BioNERO) with:

``` r
BiocManager::install("almeidasilvaf/BioNERO")
```

## Citation

Below is the citation output from using `citation('BioNERO')` in R.
Please run this yourself to check for any updates on how to cite
**BioNERO**.

``` r
print(citation('BioNERO'), bibtex = TRUE)
# 
# To cite BioNERO in publications use:
# 
#   Almeida-Silva, F., Venancio, T.M. BioNERO: an all-in-one
#   R/Bioconductor package for comprehensive and easy biological network
#   reconstruction. Funct Integr Genomics 22, 131-136 (2022).
#   https://doi.org/10.1007/s10142-021-00821-9
# 
# A BibTeX entry for LaTeX users is
# 
#   @Article{,
#     title = {BioNERO: an all-in-one R/Bioconductor package for comprehensive and easy biological network reconstruction},
#     author = {Fabricio Almeida-Silva and Thiago M. Venancio},
#     journal = {Functional And Integrative Genomics},
#     year = {2022},
#     volume = {22},
#     number = {1},
#     pages = {131-136},
#     url = {https://link.springer.com/article/10.1007/s10142-021-00821-9},
#     doi = {10.1007/s10142-021-00821-9},
#   }
```

Please note that the `BioNERO` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.
