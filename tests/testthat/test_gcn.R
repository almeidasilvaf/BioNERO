context("Test GCN inference and analysis")
library(BioNERO)

#----Simulate expression and correlation matrices----
exp <- t(matrix(rnorm(10000), ncol=1000, nrow=200))
rownames(exp) <- paste0("Gene", 1:nrow(exp))
colnames(exp) <- paste0("Sample", 1:ncol(exp))
cormat <- cor(t(exp))

#----Start tests----
test_that("get_edge_list() generates an edge list as a 3-column data frame", {
    genes <- paste0("Gene", 1:200)
    net <- list(correlation_matrix = cormat)
    n <- ncol(exp)
    edges_nofilter <- get_edge_list(net, genes=genes, filter = FALSE)
    expect_equal(ncol(edges_nofilter), 3)
})
