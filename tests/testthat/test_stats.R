
#----Load data------------------------------------------------------------------
set.seed(12)

# Simulate adjacency matrices
## GCN
exp <- matrix(
    rnorm(1000, mean = 10, sd = 2), nrow = 100, ncol = 100, byrow = TRUE
)
adj <- abs(cor(t(exp)))
rownames(adj) <- paste0("Gene", 1:100)
colnames(adj) <- paste0("Gene", 1:100)

## GRN and PPI
adj_un <- adj
adj_un[adj_un > 0.4] <- 1
adj_un[adj_un <= 0.4] <- 0
diag(adj_un) <- 0

#----Start tests----------------------------------------------------------------
test_that("net_stats() returns a list of network statistics", {

    stats_gcn <- net_stats(adj, net_type = "gcn")
    stats_ppi <- net_stats(adj_un, net_type = "ppi")
    stats_grn <- net_stats(adj_un, net_type = "grn")
    stats_gcn_additional <- net_stats(adj, net_type = "gcn", TRUE)

    expect_error(net_stats(adj, net_type = "wrong"))

    expect_equal(class(stats_gcn), "list")
    expect_equal(class(stats_ppi), "list")
    expect_equal(class(stats_grn), "list")
    expect_equal(class(stats_gcn_additional), "list")

    expect_equal(length(stats_gcn), 8)
    expect_equal(length(stats_ppi), 7)
    expect_equal(length(stats_grn), 6)
    expect_equal(length(stats_gcn_additional), 10)
})
