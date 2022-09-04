
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
    Class = paste0("Class", sample(1:5, size = ncol(exp), replace = TRUE))
)

## Gene metadata
row_metadata <- data.frame(
    row.names = rownames(exp),
    Pathway = paste0("Pathway", sample(1:5, size = nrow(exp), replace = TRUE))
)

## Edge list
edges <- igraph::get.edgelist(igraph::barabasi.game(n = 50, directed = FALSE))

#----Start tests----------------------------------------------------------------
test_that("sample_cols_heatmap() returns a list of processed data", {

    s1 <- sample_cols_heatmap(col_metadata, exp)
    s2 <- sample_cols_heatmap(
        cbind(col_metadata, col_metadata), exp
    )
    s3 <- sample_cols_heatmap(NULL, exp)

    # Do 3+ columns return error?
    expect_error(
        sample_cols_heatmap(
            cbind(col_metadata, col_metadata, col_metadata), exp
        )
    )

    expect_equal(class(s1), "list")
    expect_equal(class(s2), "list")
    expect_equal(class(s3), "list")

    expect_equal(length(s1), 3)
    expect_equal(length(s2), 3)
    expect_equal(length(s3), 3)

})


test_that("gene_cols_heatmap() returns a list of processed data", {

    s1 <- sample_cols_heatmap(col_metadata, exp)

    g1 <- gene_cols_heatmap(row_metadata, exp, row_metadata)
    g2 <- gene_cols_heatmap(row_metadata, exp, s1)

    expect_equal(class(g1), "list")
    expect_equal(class(g2), "list")

    expect_equal(length(g1), 3)
    expect_equal(length(g2), 3)

})


test_that("calculate_cor_adj() returns a list of cor matrix and adj matrix", {

    l1 <- calculate_cor_adj(cor_method = "spearman", exp, 8, "unsigned")
    l2 <- calculate_cor_adj(cor_method = "biweight", exp, 8, "unsigned")

    expect_error(
        calculate_cor_adj(cor_method = "wrong", exp, 8, "unsigned")
    )

    expect_equal(length(l1), 2)
    expect_equal(length(l2), 2)

    expect_equal(class(l1), "list")
    expect_equal(class(l2), "list")

})


test_that("get_TOMtype() returns a character of TOM type", {
    t <- get_TOMtype("unsigned")

    expect_equal(t, "unsigned")
})

test_that("handle_trait_type() handles variables for trait object", {

    h <- handle_trait_type(col_metadata, continuous_trait = TRUE)

    expect_equal(class(h), "data.frame")
})


test_that("check_SFT() returns a list with SFT stats", {

    c1 <- check_SFT(edges, "gcn")
    c2 <- check_SFT(edges, "grn")
    c3 <- check_SFT(edges, "ppi")

    expect_error(check_SFT(edges, "error"))

    expect_equal(class(c1), "list")
    expect_equal(class(c2), "list")
    expect_equal(class(c3), "list")


    expect_equal(length(c1), 6)
    expect_equal(length(c2), 6)
    expect_equal(length(c3), 6)

})


test_that("handle_metadata() handles metadata", {

    h <- handle_metadata(exp, col_metadata)

    expect_equal(class(h), "data.frame")

})

