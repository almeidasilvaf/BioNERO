
#----Load data------------------------------------------------------------------
set.seed(12)

data(zma.se)
data(zma.interpro)
data(filt.se)

# Simulate data
## Expression data
exp <- matrix(
    rnorm(1000, mean = 10, sd = 2), nrow = 100, ncol = 100, byrow = TRUE
)

rownames(exp) <- paste0("Gene", seq_len(nrow(exp)))
colnames(exp) <- paste0("Sample", seq_len(ncol(exp)))

## Correlation matrix
cormat <- abs(cor(t(exp)))
rownames(cormat) <- paste0("Gene", 1:100)
colnames(cormat) <- paste0("Gene", 1:100)

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

# Infer GCN to avoid repetition many test chunks
gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson")

#----Start tests----------------------------------------------------------------
test_that("SFT_fit() performs SFT fit test and returns a list with", {
    sft <- suppressWarnings(SFT_fit(se, cor_method = "pearson"))
    sft2 <- suppressWarnings(SFT_fit(se, cor_method = "biweight"))
    sft3 <- suppressWarnings(SFT_fit(se, cor_method = "spearman"))

    expect_error(
        suppressWarnings(SFT_fit(se, cor_method = "error"))
    )

    expect_equal(class(sft), "list")
    expect_equal(length(sft), 2)
    expect_true("patchwork" %in% class(sft$plot))
    expect_true("ggplot" %in% class(sft$plot))
    expect_equal(class(sft$power), "numeric")

    expect_equal(class(sft2), "list")
    expect_equal(class(sft3), "list")

    expect_equal(length(sft2), 2)
    expect_equal(length(sft3), 2)
})


test_that("exp2gcn() returns a list and get_hubs_gcn() finds hubs", {

    # exp2gcn()
    gcn <- exp2gcn(
        filt.se[1:100, ], SFTpower = 18, cor_method = "pearson", verbose = TRUE
    )
    gcn2 <- exp2gcn(filt.se[1:100, ], SFTpower = 3, cor_method = "spearman")
    gcn3 <- exp2gcn(filt.se[1:100, ], SFTpower = 18, cor_method = "biweight")

    expect_error(
        exp2gcn(filt.se, cor_method = "pearson")
    )
    expect_error(
        exp2gcn(filt.se[1:100, ], SFTpower = 3, cor_method = "error")
    )

    expect_equal(class(gcn), "list")
    expect_equal(length(gcn), 7)

    expect_equal(class(gcn2), "list")
    expect_equal(class(gcn3), "list")

    # get_hubs_gcn()
    hubs <- get_hubs_gcn(filt.se[1:100, ], gcn)
    hubs2 <- get_hubs_gcn(filt.se[1:100, ], gcn2)
    hubs3 <- get_hubs_gcn(filt.se[1:100, ], gcn3)


    expect_equal(ncol(hubs), 3)
    expect_equal(class(hubs), "data.frame")
    expect_equal(colnames(hubs), c("Gene", "Module", "kWithin"))

    expect_equal(class(hubs2), "data.frame")
    expect_equal(class(hubs3), "data.frame")

    expect_equal(ncol(hubs2), 3)
    expect_equal(ncol(hubs3), 3)
})

test_that("plot_eigengene_network() plots eigengene networks", {

    p <- plot_eigengene_network(gcn)
    expect_equal(attr(class(p), "package"), "ComplexHeatmap")

})

test_that("plot_dendro_and_colors() plots dendro and colors", {

    p <- plot_dendro_and_colors(gcn)
    expect_true("ggplot" %in% class(p))

})

test_that("module_stability() recomputes network with n resamplings", {

    p <- module_stability(exp = filt.se[1:50, ], net = gcn, nRuns = 1)

    expect_true("ggplot" %in% class(p))
})

test_that("module_trait_cor() and plot_module_trait_cor() work", {

    # module_trait_cor()
    mod_trait <- module_trait_cor(filt.se, MEs = gcn$MEs)

    expect_equal(class(mod_trait), "data.frame")
    expect_equal(ncol(mod_trait), 5)
    expect_equal(names(mod_trait), c("ME", "trait", "cor", "pvalue", "group"))


})


test_that("gene_significance() returns a data frame of gene-trait cors", {
    gs <- gene_significance(filt.se)
    gs2 <- gene_significance(
        filt.se, genes = rownames(filt.se)[1:5], use_abs = FALSE
    )

    p <- plot_gene_significance(gs)

    expect_equal(class(gs), "data.frame")
    expect_equal(ncol(gs), 5)

    expect_equal(class(gs2), "data.frame")
    expect_equal(ncol(gs2), 5)

    expect_equal(attr(class(p), "package"), "ComplexHeatmap")
})


test_that("enrichment_analysis() performs overrepresentation analysis", {

    # Simulated annotation with 1 class
    annotation <- data.frame(
        Gene = rownames(row_metadata),
        Class = row_metadata$Pathway
    )
    # Simulated annotation with 2 classes
    annotation2 <- data.frame(
        Gene = rownames(row_metadata),
        Class1 = row_metadata$Pathway,
        Class2 = row_metadata$Pathway
    )

    genes <- annotation$Gene[1:10]
    background_genes <- annotation$Gene

    # ora()
    genesets <- split(annotation, annotation$Class)
    pe <- ora(genes, genesets, background_genes)


    # enrichment_analysis()
    e1 <- enrichment_analysis(genes, background_genes, annotation, p = 1)
    e2 <- enrichment_analysis(
        genes, background_genes, annotation, column = 2, p = 0.001
    )
    e3 <- enrichment_analysis(
        genes, background_genes, annotation2,
        column = c("Class1", "Class2"), p = 1
    )

    e4 <- enrichment_analysis(
        genes, background_genes, annotation2,
        column = c("Class1", "Class2"), p = 0.0001
    )

    expect_equal(class(pe), "data.frame")

    expect_equal(class(e1), "data.frame")
    expect_equal(ncol(e1), 6)

    expect_true(nrow(e2) < 1)

    expect_equal(class(e3), "data.frame")
    expect_equal(ncol(e3), 6)

    expect_true(nrow(e4) < 1)
})


test_that("module_enrichment() performs ORA for all modules", {

    # Filter GCN to increase speed
    gcn2 <- gcn
    gcn2$genes_and_modules <- gcn2$genes_and_modules[
        gcn2$genes_and_modules$Modules == "black",
    ]

    background <- rownames(filt.se)

    me1 <- module_enrichment(
        net = gcn2, background_genes = background,
        annotation = zma.interpro, p = 1
    )
    me2 <- module_enrichment(
        net = gcn2, background_genes = background,
        annotation = zma.interpro, p = 0.000000001
    )

    expect_equal(class(me1), "data.frame")
    expect_true(is.null(me2))
})


test_that("get_neighbors() returns a list of neighbors for each gene", {

    genes <- rownames(filt.se)[1:10]
    gcn2 <- gcn
    gcn2$params$net_type <- "unsigned"

    neighbors <- get_neighbors(genes, gcn)
    neighbors2 <- get_neighbors(genes, gcn2)

    expect_equal(class(neighbors), "list")
    expect_equal(length(neighbors), 10)

    expect_equal(class(neighbors2), "list")
})


test_that("get_edge_list() generates an edge list as a 3-column data frame", {

    genes <- paste0("Gene", 1:100)
    net <- list(correlation_matrix = cormat)
    n <- ncol(exp)

    edges_nofilter <- get_edge_list(net, genes = genes, filter = FALSE)
    edges_filter_optimalsft <- get_edge_list(
        net, genes = genes, filter = TRUE, r_optimal_test = c(0,1, by = 0.02)
    )
    edges_filter_pval <- get_edge_list(
        net, genes = genes, filter = TRUE, method = "pvalue", nSamples = n
    )
    edges_filter_Z <- get_edge_list(
        net, genes = genes, filter = TRUE, method = "Zscore"
    )
    edges_filter_mincor <- get_edge_list(
        net, genes = genes, filter = TRUE, method = "min_cor", rcutoff = 0.1
    )

    expect_error(
        get_edge_list(
            net, genes = genes, filter = TRUE, method = "pvalue"
        )
    )

    expect_error(
        get_edge_list(
            net, genes = genes, filter = TRUE, method = "error", nSamples = n
        )
    )

    expect_equal(ncol(edges_nofilter), 3)
    expect_equal(ncol(edges_filter_optimalsft), 3)
    expect_equal(ncol(edges_filter_pval), 3)
    expect_equal(ncol(edges_filter_Z), 3)
    expect_equal(ncol(edges_filter_mincor), 3)
})
