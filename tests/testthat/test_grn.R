
#----Simulate expression matrix----
set.seed(123)
exp <- matrix(rnorm(10000), ncol=200, nrow=1000)
rownames(exp) <- paste0("Gene", 1:nrow(exp))
colnames(exp) <- paste0("Sample", 1:ncol(exp))

#----Infer GRNs----
clr <- grn_clr(exp)
aracne <- grn_aracne(exp)
genie3 <- grn_genie3(exp=exp, nTrees=5)
grn_list <- grn_combined(exp, nTrees=5)
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
    expect_equal(all.equal(nrow1, nrow2, nrow3), TRUE)
    expect_equal(class(grn_list), "list")
    expect_equal(length(grn_list), 3)
})


test_that("grn_average_rank() ranks GRN weights and calculate the average across methods", {
    row_test <- sample(seq_len(nrow(ranked_grn)), size=1)
    expect_equal(ncol(ranked_grn), 3)
    expect_equal(ranked_grn[row_test, 3] < ranked_grn[row_test+1, 3], TRUE)
})


test_that("check_sft() checks SFT fit for different network types", {
    expect_message(check_sft(clr, net_type = "grn"), "Your graph fits the scale-free topology. P-value:0.905436442816135")
})


test_that("grn_filter() calculates the best SFT fit for n number of top edges", {
    filtered_edges <- grn_filter(ranked_grn, nsplit=5)
    expect_equal(ncol(filtered_edges), 2)
    expect_lte(nrow(filtered_edges), nrow(ranked_grn))
    expect_message(grn_filter(ranked_grn, nsplit=5), "The top number of edges that best fits the scale-free topology is 90900")
})









