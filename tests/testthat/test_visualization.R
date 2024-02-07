
#----Load data------------------------------------------------------------------
data(filt.se)
tfs <- sample(rownames(filt.se), size=50, replace=FALSE)

## Infer networks to be used
set.seed(123)
ppi_edges <- igraph::as_edgelist(igraph::sample_pa(n=50))
gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson")
gcn_edges <- get_edge_list(
    gcn, module = "brown", filter = TRUE, method = "min_cor"
)


#----Start tests----------------------------------------------------------------
test_that("detect_communities() returns a data frame of communities", {
    com <- detect_communities(ppi_edges, directed=TRUE)
    com2 <- detect_communities(ppi_edges, directed=FALSE)
    expect_equal(class(com), "data.frame")
    expect_equal(class(com2), "data.frame")
})


test_that("plot_ppi() plots a PPI network", {

    color_by <- data.frame(
        Node = 1:50,
        Class = paste0("Class", sample(c(1:50), size = 50, replace = TRUE))
    )

    p <- plot_ppi(ppi_edges, add_color_legend = FALSE)
    p2 <- plot_ppi(ppi_edges, color_by = color_by)
    p3 <- plot_ppi(ppi_edges, add_color_legend = FALSE, show_labels = "all")
    p4 <- plot_ppi(ppi_edges, add_color_legend = FALSE, show_labels = "allhubs")
    p5 <- plot_ppi(ppi_edges, add_color_legend = FALSE, show_labels = "none")

    expect_error(
        plot_ppi(ppi_edges, add_color_legend = FALSE, show_labels = "any")
    )

    pint <- plot_ppi(ppi_edges, interactive = TRUE)

    expect_equal(class(p), c("gg", "ggplot"))
    expect_equal(class(p2), c("gg", "ggplot"))
    expect_equal(class(p3), c("gg", "ggplot"))
    expect_equal(class(p4), c("gg", "ggplot"))
    expect_equal(class(p5), c("gg", "ggplot"))

    expect_true("forceNetwork" %in% class(pint))

})


test_that("plot_grn() plots a GRN", {

    p <- plot_grn(ppi_edges, ranked = FALSE)
    p3 <- plot_grn(ppi_edges, show_labels = "all")
    p4 <- plot_grn(ppi_edges, show_labels = "allhubs")
    p5 <- plot_grn(ppi_edges, show_labels = "none")

    expect_error(
        plot_ppi(ppi_edges, show_labels = "any")
    )

    pint <- plot_grn(ppi_edges, interactive = TRUE)


    expect_equal(class(p), c("gg", "ggplot"))
    expect_equal(class(p3), c("gg", "ggplot"))
    expect_equal(class(p4), c("gg", "ggplot"))
    expect_equal(class(p5), c("gg", "ggplot"))

    expect_true("forceNetwork" %in% class(pint))

})


test_that("plot_gcn() plots a GCN", {

    color_by <- data.frame(
        Node = unique(c(gcn_edges$Gene1, gcn_edges$Gene2)),
        Class = paste0("Class", sample(c(1:50), size = 32, replace = TRUE))
    )
    hubs <- data.frame(
        Gene = c("ZeamMp044", "ZeamMp092", "ZeamMp108"),
        Degree = rep(10, 3)
    )


    p <- plot_gcn(gcn_edges, gcn, hubs = hubs)
    p2 <- plot_gcn(gcn_edges, gcn, hubs = hubs, color_by = color_by)
    p3 <- plot_gcn(gcn_edges, gcn, hubs = hubs, show_labels = "all")
    p4 <- plot_gcn(gcn_edges, gcn, hubs = hubs, show_labels = "allhubs")
    p5 <- plot_gcn(gcn_edges, gcn, hubs = hubs, show_labels = "none")

    expect_error(
        plot_gcn(gcn_edges)
    )

    expect_error(
        plot_gcn(gcn_edges, gcn, hubs = hubs, show_labels = "wrong")
    )

    pint <- plot_gcn(gcn_edges, gcn, hubs = hubs, show_labels = "allhubs",
                     interactive = TRUE)

    expect_equal(class(p), c("gg", "ggplot"))
    expect_equal(class(p2), c("gg", "ggplot"))
    expect_equal(class(p3), c("gg", "ggplot"))
    expect_equal(class(p4), c("gg", "ggplot"))
    expect_equal(class(p5), c("gg", "ggplot"))

    expect_true("forceNetwork" %in% class(pint))
})
