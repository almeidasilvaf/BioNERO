context("Expression data processing")
library(BioNERO)

test_that("dfs2one reads in multiple tables and binds them into a data frame", {
    genes <- paste0(rep("Gene", 100), 1:100)
    samples1 <- paste0(rep("Sample", 30), 1:30)
    samples2 <- paste0(rep("Sample", 30), 31:60)
    # Simulate an expression data frame of 100 genes and 30 samples
    exp1 <- as.data.frame(matrix(rnorm(100*30),nrow=100,ncol=30))
    colnames(exp1) <- samples1
    exp1 <- cbind(genes, exp1)

    # Simulate an expression data frame of 100 genes and 30 different samples
    exp2 <- as.data.frame(matrix(rnorm(100*30),nrow=100,ncol=30))
    colnames(exp2) <- samples2
    exp2 <- cbind(genes, exp2)

    write.table(exp1, file="exp1.tsv", quote=FALSE, sep="\t")
    write.table(exp2, file="exp2.tsv", quote=FALSE, sep="\t")

    expect_equal(nrow(dfs2one(mypath="./")), 100)
    expect_equal(ncol(dfs2one(mypath="./")), 60)

    unlink("exp1.tsv")
    unlink("exp2.tsv")
})

