
#----Load data------------------------------------------------------------------
set.seed(12)

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
    Class = paste0("Class", sample(1:5, size = ncol(exp), replace = TRUE)),
    Weight = stats::rnorm(ncol(exp), 50, 15)
)

## Gene metadata
row_metadata <- data.frame(
    row.names = rownames(exp),
    Pathway = paste0("Pathway", sample(1:5, size = nrow(exp), replace = TRUE))
)

## Simulate a SummarizedExperiment object
se <- SummarizedExperiment::SummarizedExperiment(
    assays = exp, colData = col_metadata, rowData = row_metadata
)

## Edge list
edges <- igraph::as_edgelist(igraph::sample_pa(n = 50))

#----Start tests----------------------------------------------------------------
test_that("metadata2colors() returns a list of metadata and named vectors", {

    cols <- metadata2colors(col_metadata)

    expect_equal(names(cols), c("metadata", "colors"))
    expect_equal(length(cols), 2)

    metadata_4columns <- cbind(
        col_metadata, col_metadata, col_metadata, col_metadata
    )

    expect_error(metadata2colors(metadata_4columns))
})


test_that("heatmap_attributes() returns heatmap parameters as a list", {

    p <- heatmap_attributes(exp)

    expect_equal(length(p), 4)
    expect_equal(names(p), c("pal", "mat", "title", "name"))

    expect_error(heatmap_attributes(exp, heatmap_type = "wrong"))
})


test_that("se2metadata() returns a list of row and coldata", {

    m <- se2metadata(se)

    expect_equal(length(m), 2)
    expect_equal(names(m), c("rowdata", "coldata"))
})

test_that("get_model_matrix() returns a model matrix for module-trait cor", {

    mat_char <- get_model_matrix(col_metadata, 1)
    mat_num <- get_model_matrix(col_metadata, 2)

    expect_equal(ncol(mat_char), 5)
    expect_equal(ncol(mat_num), 1)
})


test_that("exp2cor() and cor2adj() return symmetric matrices", {

    exp <- matrix(
        rnorm(100 * 50, mean = 10, sd = 2),
        nrow = 100,
        dimnames = list(
            paste0("gene", seq_len(100)),
            paste0("sample", seq_len(50))
        )
    )
    cor_mat <- exp2cor(exp)
    cor_mat2 <- exp2cor(exp, cor_method = "spearman")
    cor_mat3 <- exp2cor(exp, cor_method = "biweight")
    adj <- cor2adj(cor_mat, beta = 9)


    expect_error(cor2adj(cor_mat, beta = 9, net_type = "error"))
    expect_error(exp2cor(exp, cor_method = "error"))

    expect_true(is(cor_mat, "matrix"))
    expect_true(is(cor_mat2, "matrix"))
    expect_true(is(cor_mat3, "matrix"))
    expect_true(is(adj, "matrix"))

})


test_that("guess_tom() returns a character with TOM type", {

    t <- guess_tom("unsigned")

    expect_equal(t, "unsigned")
})


test_that("check_SFT() returns a list with SFT stats", {

    c1 <- check_SFT(edges, "gcn")
    c2 <- check_SFT(edges, "grn")
    c3 <- check_SFT(edges, "ppi")

    expect_error(check_SFT(edges, "error"))

    expect_equal(class(c1), "list")
    expect_equal(class(c2), "list")
    expect_equal(class(c3), "list")

})


test_that("handle_metadata() handles metadata", {

    h <- handle_metadata(exp, col_metadata)

    expect_equal(class(h), "data.frame")

})

