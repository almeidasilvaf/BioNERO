
#----Simulate expression and correlation matrices----
set.seed(1)
exp <- t(matrix(rnorm(10000), ncol=1000, nrow=200))
rownames(exp) <- paste0("Gene", 1:nrow(exp))
colnames(exp) <- paste0("Sample", 1:ncol(exp))
cormat <- cor(t(exp))

# Load data
data(zma.se)
data(zma.interpro)
data(filt.se)

# Infer GCN to avoid repetition many test chunks
gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson",
               reportPDF = FALSE)

#----Start tests----
test_that("SFT_fit() performs SFT fit test and returns a list with", {
    sft <- SFT_fit(filt.se, cor_method="pearson")
    expect_equal(class(sft), "list")
    expect_equal(length(sft), 2)
    expect_true(all.equal(class(sft$plot), c("gg", "ggplot", "ggarrange")))
    expect_equal(class(sft$power), "numeric")
})


test_that("exp2gcn() infers GCN and returns as a list", {
    gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson",
                   reportPDF = FALSE)
    expect_equal(class(gcn), "list")
    expect_equal(length(gcn), 7)
})


test_that("module_stability() recomputes network with n resamplings", {
    module_stability(exp = filt.se, net = gcn, nRuns = 1)
    output_file <- paste0(Sys.Date(), "_module_stability.pdf")
    expect_true(file.exists(output_file))
    unlink(output_file)
})


test_that("module_trait_cor() returns a heatmap of mod-trait correlations", {
    mod_trait <- module_trait_cor(filt.se, MEs=gcn$MEs)
    expect_equal(class(mod_trait), "data.frame")
    expect_equal(ncol(mod_trait), 4)
    expect_equal(class(mod_trait$ME), "character")
    expect_equal(class(mod_trait$trait), "character")
})


test_that("gene_significance() returns a list of GS matrices", {
    gs <- gene_significance(filt.se)
    expect_equal(class(gs), "list")
    expect_equal(length(gs), 2)
})


test_that("get_hubs_gcn() identified GCN hubs", {
    hubs <- get_hubs_gcn(filt.se, gcn)
    expect_equal(ncol(hubs), 3)
    expect_equal(class(hubs), "data.frame")
    expect_equal(colnames(hubs), c("Gene", "Module", "kWithin"))
})


test_that("enrichment_analysis() performs overrepresentation analysis", {
    BiocParallel::register(BiocParallel::SerialParam())
    genes <- rownames(filt.se)[1:50]
    background_genes <- rownames(filt.se)
    annotation <- zma.interpro
    enrich <- enrichment_analysis(genes, background_genes,
                                  annotation, p = 1)
    expect_equal(class(enrich), "data.frame")
    expect_equal(ncol(enrich), 6)
})

# Test passed
# The code was kept here so we can return to this test if necessary.
# As the test took a long time to run, we decided to comment it.
#
# test_that("module_enrichment() performs ORA for all modules", {
#     BiocParallel::register(BiocParallel::SerialParam())
#     background <- rownames(filt.se)
#     mod_enrich <- module_enrichment(net = gcn,
#                                     background_genes = background,
#                                     annotation = zma.interpro, p=1)
#     expect_equal(class(mod_enrich), "data.frame")
# })


test_that("get_neighbors() returns a list of neighbors for each gene", {
    genes <- rownames(filt.se)[1:10]
    neighbors <- get_neighbors(genes, gcn)
    expect_equal(class(neighbors), "list")
    expect_equal(length(neighbors), 10)
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
