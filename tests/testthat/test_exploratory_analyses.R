
#----Load data----------------------------------------------------------------
set.seed(12)
data(filt.se)

# Simulate data
## Expression data
exp <- matrix(
    rnorm(1000, mean = 10, sd = 2), nrow = 100, ncol = 100, byrow = TRUE
)
rownames(exp) <- paste0("Gene", seq_len(nrow(exp)))
colnames(exp) <- paste0("Sample", seq_len(ncol(exp)))

## Sample metadata
col_metadata <- data.frame(
    row.names = colnames(exp),
    Class = paste0("Class", sample(1:5, size = ncol(exp), replace = TRUE))
)

## Gene metadata
row_metadata <- data.frame(
    row.names = rownames(exp),
    Pathway = paste0("Pathway", sample(1:5, size = nrow(exp), replace = TRUE))
)

## SummarizedExperiment object
se <- SummarizedExperiment::SummarizedExperiment(
    exp, colData = col_metadata
)

# Get GCN
gcn <- suppressWarnings(
    exp2gcn(filt.se, SFTpower = 16, cor_method = "pearson")
)

#----Start tests----------------------------------------------------------------
test_that("plot_heatmap() correctly handles row and col annotation", {

    p <- plot_heatmap(
        exp, col_metadata, row_metadata, cluster_cols = FALSE,
        cluster_rows = FALSE, type = "expr", log_trans = TRUE
    )

    p2 <- plot_heatmap(
        se, col_metadata, row_metadata, cluster_cols = FALSE,
        cluster_rows = FALSE, type = "samplecor", log_trans = FALSE,
        palette = "Greens"
    )

    expect_error(
        plot_heatmap(
            exp, col_metadata, row_metadata, cluster_cols = FALSE,
            cluster_rows = FALSE, type = "error", log_trans = TRUE
        )
    )


    expect_equal(attr(class(p), "package"), "ComplexHeatmap")
    expect_equal(attr(class(p2), "package"), "ComplexHeatmap")

})


test_that("plot_PCA() performs PCA and plots it", {

    p1 <- plot_PCA(exp, col_metadata, log_trans = TRUE)
    p2 <- plot_PCA(filt.se, log_trans = FALSE, PCs = c(1,3))

    expect_error(
        plot_PCA(exp, cbind(col_metadata, col_metadata, col_metadata))
    )

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

    p1 <- plot_expression_profile(
        genes = genes, exp = zma.se, plot_module = FALSE
    )

    p2 <- plot_expression_profile(
        exp = filt.se, plot_module = TRUE, net = gcn, modulename = "black"
    )

    expect_true("ggplot" %in% class(p1))
    expect_true("ggplot" %in% class(p2))
})


test_that("plot_ngenes_per_module() returns a barplot of genes per module", {
    p <- plot_ngenes_per_module(gcn)
    expect_equal(class(p), c("gg", "ggplot"))
})


