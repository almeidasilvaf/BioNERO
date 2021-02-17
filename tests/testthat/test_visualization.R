
#----Load data----
data(filt.se)
tfs <- sample(rownames(filt.se), size=50, replace=FALSE)

#---Infer networks to be used----
set.seed(123)
ppi_edges <- igraph::get.edgelist(igraph::barabasi.game(n=50, directed=FALSE))
grn_edges <- grn_clr(filt.se, regulators = tfs)
gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson",
                              reportPDF = FALSE)
gcn_edges <- get_edge_list(gcn, module="brown", filter=TRUE, method="min_cor")


#----Start tests----
test_that("detect_communities() returns a data frame of communities", {
    com <- detect_communities(grn_edges, directed=TRUE)
    com2 <- detect_communities(ppi_edges, directed=FALSE)
    expect_equal(class(com), "data.frame")
    expect_equal(class(com2), "data.frame")
})


test_that("plot_ppi() plots a PPI network", {
    p <- plot_ppi(ppi_edges, add_color_legend = FALSE)
    expect_equal(class(p), c("gg", "ggplot"))
})


test_that("plot_grn() plots a GRN", {
    grn_final <- exp2grn(filt.se, regulators = tfs, nTrees=5)
    p <- plot_grn(grn_edges, ranked=FALSE)
    expect_equal(class(p), c("gg", "ggplot"))
})


test_that("plot_gcn() plots a GCN", {
    gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson",
                   reportPDF = FALSE)
    gcn_edges <- get_edge_list(gcn, module="brown", filter=TRUE,
                               method="min_cor")
    hubs <- get_hubs_gcn(filt.se, gcn)
    p <- plot_gcn(gcn_edges, gcn, hubs = hubs)
    expect_equal(class(p), c("gg", "ggplot"))
})
