
#----Create data for testing----
set.seed(12)
data("zma.se")
filt.zma <- filter_by_variance(zma.se, n=500)
zma.set1 <- filt.zma[, sample(colnames(filt.zma), size=20, replace=FALSE)]
zma.set2 <- filt.zma[, sample(colnames(filt.zma), size=20, replace=FALSE)]
list.sets <- list(zma.set1, zma.set2)

#----Start tests----

test_that("consensus_SFT_fit() calculates best SFT fit for all networks", {
    cons_sft <- consensus_SFT_fit(list.sets,
                                  setLabels = c("Maize 1", "Maize 2"),
                                  cor_method = "pearson")
    expect_equal(class(cons_sft[[1]]), "numeric")
    expect_equal(class(cons_sft[[2]]), c("gg", "ggplot", "ggarrange"))
    expect_equal(length(cons_sft[[1]]), length(list.sets))
})

test_that("consensus_modules() identified consensus modules across sets", {
    cons_mod <- consensus_modules(list.sets, power = c(11, 13), cor_method = "pearson")
    expect_equal(class(cons_mod), "list")
    expect_equal(length(cons_mod), 6)
})

test_that("plot_dendro_and_cons_colors() plots dendro and colors", {
    cons_mod <- consensus_modules(list.sets, power = c(11, 13), cor_method = "pearson")
    p <- plot_dendro_and_cons_colors(cons_mod)
    expect_equal(class(p), "list")
})

test_that("consensus_trait_cor() correlates consensus mods to traits", {
    consensus <- consensus_modules(list.sets, power = c(11, 13), cor_method = "pearson")
    consensus_trait <- consensus_trait_cor(consensus, cor_method = "pearson")
    expect_equal(class(consensus_trait), "data.frame")
    expect_equal(names(consensus_trait), c("ME", "trait", "cor", "pvalue"))
})

