context("GRN inference")
library(BioNERO)

#----Simulate expression matrix----
exp <- t(matrix(rnorm(10000), ncol=1000, nrow=200))
rownames(exp) <- paste0("Gene", 1:nrow(exp))
colnames(exp) <- paste0("Sample", 1:ncol(exp))

#----Start tests----
test_that("cormat_to_edgelist() converts correlation matrix to edge list", {
    cor_mat <- cor(t(exp))
    edgelist <- cormat_to_edgelist(cor_mat)
    expect_equal(ncol(edgelist), 3)
    expect_equal(colnames(edgelist), c("Node1", "Node2", "Weight"))
})


test_that("grn_clr() produces an edge list", {
    grn <- grn_clr(exp)
    expect_equal(colnames(grn), c("Node1", "Node2", "Weight"))
    expect_equal(is.character(grn[,1]), TRUE)
    expect_equal(is.character(grn[,2]), TRUE)
    expect_equal(is.numeric(grn[,3]), TRUE)
})


test_that("grn_aracne() produces an edge list", {
    grn <- grn_aracne(exp)
    expect_equal(ncol(grn), 3)
    expect_equal(colnames(grn), c("Node1", "Node2", "Weight"))
    expect_equal(is.character(grn[,1]), TRUE)
    expect_equal(is.character(grn[,2]), TRUE)
    expect_equal(is.numeric(grn[,3]), TRUE)
})


test_that("grn_genie3() produces an edge list", {
    grn <- grn_genie3(exp=exp, nTrees=100)
    expect_equal(ncol(grn), 3)
    expect_equal(is.character(grn[,1]), TRUE)
    expect_equal(is.character(grn[,2]), TRUE)
    expect_equal(is.numeric(grn[,3]), TRUE)
})
















