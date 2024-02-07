
#---Load data----
set.seed(1)
data(filt.se)
tfs <- sample(rownames(filt.se), size=50, replace=FALSE)

#----Infer GRNs----
clr <- grn_infer(filt.se, method = "clr", regulators = tfs)
aracne <- grn_infer(filt.se, method = "aracne", regulators = tfs)
genie3 <- grn_infer(filt.se, method = "genie3", regulators = tfs, nTrees = 1)
grn_list <- grn_combined(filt.se, regulators = tfs, nTrees = 1)
ranked_grn <- grn_average_rank(grn_list)

#----Start tests----
test_that("cormat_to_edgelist() converts correlation matrix to edge list", {
    cor_mat <- cor(t(SummarizedExperiment::assay(filt.se)))
    edgelist <- cormat_to_edgelist(cor_mat)
    expect_equal(ncol(edgelist), 3)
    expect_equal(colnames(edgelist), c("Node1", "Node2", "Weight"))
})

test_that("grn_infer() produces an edge list", {

    expect_equal(colnames(clr), c("Node1", "Node2", "Weight"))
    expect_equal(is.character(clr[,1]), TRUE)
    expect_equal(is.character(clr[,2]), TRUE)
    expect_equal(is.numeric(clr[,3]), TRUE)
    expect_equal(ncol(aracne), 3)
    expect_equal(colnames(aracne), c("Node1", "Node2", "Weight"))
    expect_equal(is.character(aracne[,1]), TRUE)
    expect_equal(is.character(aracne[,2]), TRUE)
    expect_equal(is.numeric(aracne[,3]), TRUE)
    expect_equal(ncol(genie3), 3)
    expect_equal(is.character(genie3[,1]), TRUE)
    expect_equal(is.character(genie3[,2]), TRUE)
    expect_equal(is.numeric(genie3[,3]), TRUE)

    expect_error(grn_infer(filt.se, method = "random", regulators = tfs))
    expect_error(grn_infer(filt.se, method = "clr"))
})

test_that("grn_combined() produces a list of edge lists", {
    nrow1 <- nrow(grn_list[[1]])
    nrow2 <- nrow(grn_list[[2]])
    nrow3 <- nrow(grn_list[[3]])
    expect_equal(class(grn_list), "list")
    expect_equal(length(grn_list), 3)
})

test_that("all inference algorithms return regulators in the first column
          and target in the second column", {
    edge_structure <- function(grn, tfs) {
        x <- nrow(grn[grn[,2] %in% tfs, ])
        y <- nrow(grn[grn[,1] %in% tfs, ])
        result <- c(x,y)
        return(result)
    }
    expect_equal(edge_structure(genie3, tfs), c(0, nrow(genie3)))
    expect_equal(edge_structure(clr, tfs), c(0, nrow(clr)))
    expect_equal(edge_structure(aracne, tfs), c(0, nrow(aracne)))
})


test_that("grn_average_rank() ranks GRN weights and calculate the average across methods", {
    row_test <- sample(seq_len(nrow(ranked_grn)), size=1)
    expect_equal(ncol(ranked_grn), 3)
    expect_equal(ranked_grn[row_test, 3] <= ranked_grn[row_test+1, 3], TRUE)
})


test_that("check_SFT() checks SFT fit for different network types", {
    check <- check_SFT(clr, net_type = "grn")
    expect_equal(class(check), "list")
})


test_that("grn_filter() calculates the best SFT fit for n number of top edges", {
    filtered_edges <- grn_filter(ranked_grn, nsplit = 2)

    expect_equal(ncol(filtered_edges), 2)
    expect_lte(nrow(filtered_edges), nrow(ranked_grn))
})


test_that("get_hubs_grn() returns a data frame with hub genes and their out degree", {
    filtered_edges <- grn_filter(ranked_grn, nsplit = 2)
    n_genes <- length(unique(
        c(as.character(filtered_edges[,1]), as.character(filtered_edges[,2]))
    ))

    fedges <- filtered_edges
    fedges$rank <- seq_len(nrow(fedges))

    hubs <- get_hubs_grn(filtered_edges)
    hubs2 <- get_hubs_grn(filtered_edges, top_n = 1)
    hubs3 <- get_hubs_grn(fedges, top_n = 1)

    expect_equal(nrow(hubs), floor(n_genes * 0.1))
    expect_equal(nrow(hubs2), 1)
    expect_equal(nrow(hubs3), 1)
})

test_that("get_hubs_ppi() returns a list or data frame of hubs and degree", {
    ppi_edges <- igraph::as_edgelist(igraph::sample_pa(n=500))
    hubs <- get_hubs_ppi(ppi_edges, return_degree = TRUE)
    hubs2 <- get_hubs_ppi(ppi_edges, top_n = 1)

    expect_equal(class(hubs[[1]]), "data.frame")
    expect_equal(class(hubs), "list")
    expect_equal(nrow(hubs2), 1)
})


test_that("exp2grn() infers a GRN from expression", {
    grn <- exp2grn(filt.se, regulators = tfs, nTrees=2, nsplit=2)

    expect_equal(ncol(grn), 2)
    expect_equal(class(grn), "data.frame")
})







