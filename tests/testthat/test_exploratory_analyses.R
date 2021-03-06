
#----Set up data----
data(zma.se)
data(filt.se)
exp <- SummarizedExperiment::assay(zma.se)
exp <- filter_by_variance(exp, n=100)
filtered.se <- filter_by_variance(zma.se, n=100)
metadata <- as.data.frame(SummarizedExperiment::colData(zma.se))

#----Start tests----
test_that("plot_heatmap() correctly handles row and col annotation", {
    metadata2 <- metadata
    metadata2$Annot2 <- sample(c("Class1", "Class2", "Class4", "Class5",
                               "Class6"),
                               size=nrow(metadata2), replace=TRUE)
    row_annot <- data.frame(row.names=rownames(exp),
                            Class=sample(c("Pathway1", "Pathway2",
                                           "Pathway3"),
                                         size=nrow(exp), replace=TRUE))
    p <- plot_heatmap(exp, col_metadata = metadata2, row_metadata = row_annot,
                 cluster_cols=FALSE, cluster_rows=FALSE, type="expr",
                 log_trans=TRUE)
    expect_true(class(p) == "HeatmapList")
})


test_that("plot_PCA() performs PCA and plots it", {
    p1 <- plot_PCA(exp, metadata, log_trans = TRUE)
    p2 <- plot_PCA(filt.se, log_trans = TRUE)

    expect_true(all.equal(class(p1), c("gg", "ggplot")))
    expect_true(all.equal(class(p2), c("gg", "ggplot")))
})


test_that("get_HK() returns housekeeping genes in a character vector", {
    hk <- get_HK(zma.se)

    expect_equal(class(hk), "character")
    expect_true(length(hk) <= floor(nrow(zma.se) * 0.25))
})


test_that("plot_expression_profile() plots expression in a line plot", {
    genes <- rownames(filt.se)
    p1 <- plot_expression_profile(genes=genes, exp=zma.se, plot_module=FALSE)
    expect_true(all.equal(class(p1), c("gg", "ggplot")))
})


test_that("plot_ngenes_per_module() returns a barplot of genes per module", {
    gcn <- exp2gcn(filt.se, SFTpower = 16, cor_method = "pearson")
    p <- plot_ngenes_per_module(gcn)
    expect_equal(class(p), c("gg", "ggplot"))
})


