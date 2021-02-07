
#----Simulate expression and correlation matrices----
set.seed(1)
exp <- t(matrix(rnorm(10000), ncol=1000, nrow=200))
rownames(exp) <- paste0("Gene", 1:nrow(exp))
colnames(exp) <- paste0("Sample", 1:ncol(exp))
cormat <- cor(t(exp))

data(se.seed)
filt.se <- filter_by_variance(se.seed, n=500)

#----Start tests----
test_that("SFT_fit() performs SFT fit test and returns a list with", {
    sft <- SFT_fit(filt.se, cor_method="pearson")
    expect_equal(class(sft), "list")
    expect_equal(length(sft), 2)
    expect_true(all.equal(class(sft$plot), c("gg", "ggplot", "ggarrange")))
    expect_equal(class(sft$power), "numeric")
})



test_that("get_edge_list() generates an edge list as a 3-column data frame", {
    genes <- paste0("Gene", 1:200)
    net <- list(correlation_matrix = cormat)
    n <- ncol(exp)
    edges_nofilter <- get_edge_list(net, genes=genes, filter = FALSE)
    edges_filter_optimalsft <- get_edge_list(net, genes=genes, filter=TRUE, r_optimal_test = c(0,1,by=0.02))
    edges_filter_pval <- get_edge_list(net, genes=genes, filter=TRUE, method="pvalue", nSamples = n)
    edges_filter_Z <- get_edge_list(net, genes=genes, filter=TRUE, method="Zscore")
    edges_filter_mincor <- get_edge_list(net, genes=genes, filter=TRUE, method="min_cor", rcutoff = 0.1)

    expect_equal(ncol(edges_nofilter), 3)
    expect_message(
        get_edge_list(net, genes=genes, filter=TRUE, r_optimal_test = c(0,1,by=0.02)),
        "The correlation threshold that best fits the scale-free topology is 0.02"
    )
    expect_equal(ncol(edges_filter_optimalsft), 3)
    expect_equal(ncol(edges_filter_pval), 3)
    expect_equal(ncol(edges_filter_Z), 3)
    expect_equal(ncol(edges_filter_mincor), 3)
})
