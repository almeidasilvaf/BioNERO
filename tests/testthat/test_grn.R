
#----Simulate expression matrix----
set.seed(123)
exp <- matrix(rnorm(10000), ncol=200, nrow=1000)
rownames(exp) <- paste0("Gene", 1:nrow(exp))
colnames(exp) <- paste0("Sample", 1:ncol(exp))
tfs <- paste0("Gene", 1:100)

#----Infer GRNs----
clr <- grn_clr(exp, regulators = tfs)
aracne <- grn_aracne(exp, regulators = tfs)
genie3 <- grn_genie3(exp=exp, regulators=tfs, nTrees=5)
grn_list <- grn_combined(exp, regulators=tfs, nTrees=5)
ranked_grn <- grn_average_rank(grn_list)

#----Start tests----
test_that("cormat_to_edgelist() converts correlation matrix to edge list", {
    cor_mat <- cor(t(exp))
    edgelist <- cormat_to_edgelist(cor_mat)
    expect_equal(ncol(edgelist), 3)
    expect_equal(colnames(edgelist), c("Node1", "Node2", "Weight"))
})


test_that("grn_clr() produces an edge list", {
    expect_equal(colnames(clr), c("Node1", "Node2", "Weight"))
    expect_equal(is.character(clr[,1]), TRUE)
    expect_equal(is.character(clr[,2]), TRUE)
    expect_equal(is.numeric(clr[,3]), TRUE)
})


test_that("grn_aracne() produces an edge list", {
    expect_equal(ncol(aracne), 3)
    expect_equal(colnames(aracne), c("Node1", "Node2", "Weight"))
    expect_equal(is.character(aracne[,1]), TRUE)
    expect_equal(is.character(aracne[,2]), TRUE)
    expect_equal(is.numeric(aracne[,3]), TRUE)
})


test_that("grn_genie3() produces an edge list", {
    expect_equal(ncol(genie3), 3)
    expect_equal(is.character(genie3[,1]), TRUE)
    expect_equal(is.character(genie3[,2]), TRUE)
    expect_equal(is.numeric(genie3[,3]), TRUE)
})


test_that("grn_combined() produces a list of edge lists", {
    nrow1 <- nrow(grn_list[[1]])
    nrow2 <- nrow(grn_list[[2]])
    nrow3 <- nrow(grn_list[[3]])
    expect_equal(class(grn_list), "list")
    expect_equal(length(grn_list), 3)
})

test_that("all inference algorithms return regulators in the first column and target in the second column", {
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
    expect_equal(ranked_grn[row_test, 3] < ranked_grn[row_test+1, 3], TRUE)
})


test_that("check_sft() checks SFT fit for different network types", {
    expect_message(check_sft(clr, net_type = "grn"), "Your graph fits the scale-free topology. P-value:0.999964664866023")
})


test_that("grn_filter() calculates the best SFT fit for n number of top edges", {
    filtered_edges <- grn_filter(ranked_grn, nsplit=5)
    expect_equal(ncol(filtered_edges), 2)
    expect_lte(nrow(filtered_edges), nrow(ranked_grn))
})


test_that("get_hubs_grn() return a data frame with hub genes and their out degree", {
    filtered_edges <- grn_filter(ranked_grn, nsplit = 5)
    n_genes <- length(unique(c(as.character(filtered_edges[,1]), as.character(filtered_edges[,2]))))
    hubs <- get_hubs_grn(filtered_edges)
    expect_equal(nrow(hubs), floor(n_genes * 0.1))
})








