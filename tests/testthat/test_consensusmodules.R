
#----Create data for testing----------------------------------------------------
set.seed(12)
data("zma.se")
filt.zma <- filter_by_variance(zma.se, n=500)
zma.set1 <- filt.zma[, sample(colnames(filt.zma), size=20, replace=FALSE)]
zma.set2 <- filt.zma[, sample(colnames(filt.zma), size=20, replace=FALSE)]
list.sets <- list(zma.set1, zma.set2)

rlist.sets <- lapply(list.sets, function(x) return(x[1:50, ]))
rlist.sets2 <- lapply(list.sets, function(x) return(x[1:100, ]))


#----Start tests----------------------------------------------------------------
test_that("consensus_SFT_fit() calculates best SFT fit for all networks", {
    cons_sft <- consensus_SFT_fit(rlist.sets, cor_method = "pearson")
    cons_sft2 <- consensus_SFT_fit(rlist.sets, cor_method = "spearman")
    cons_sft3 <- consensus_SFT_fit(rlist.sets, cor_method = "biweight")

    expect_error(
        consensus_SFT_fit(rlist.sets, cor_method = "error")
    )

    expect_equal(class(cons_sft[[1]]), "numeric")
    expect_true("patchwork" %in% class(cons_sft[[2]]))
    expect_true("ggplot" %in% class(cons_sft[[2]]))
    expect_equal(length(cons_sft[[1]]), length(list.sets))

    expect_equal(class(cons_sft2), "list")
    expect_equal(class(cons_sft3), "list")
})

test_that("consensus_modules() and related functions work", {

    # consensus_modules()
    cm1 <- consensus_modules(
        rlist.sets2, power = c(20, 20), cor_method = "pearson", verbose = TRUE
    )
    cm2 <- consensus_modules(
        rlist.sets2, power = c(20, 20), cor_method = "spearman", verbose = TRUE
    )
    cm3 <- consensus_modules(
        rlist.sets2, power = c(15, 15), cor_method = "biweight", verbose = TRUE
    )

    expect_error(
        consensus_modules(
            rlist.sets2, power = c(20, 20), cor_method = "error"
        )

    )

    expect_equal(class(cm1), "list")
    expect_equal(length(cm1), 6)

    expect_equal(class(cm2), "list")
    expect_equal(class(cm3), "list")

    # plot_dendro_and_cons_colors()
    p <- plot_dendro_and_cons_colors(cm1)
    expect_equal(class(p), "list")

    # consensus_trait_cor()
    ct <- consensus_trait_cor(cm1, cor_method = "pearson")
    ct2 <- consensus_trait_cor(cm1, cor_method = "spearman", transpose = TRUE)

    expect_equal(class(ct), "data.frame")
    expect_equal(names(ct), c("ME", "trait", "cor", "pvalue"))

    expect_equal(class(ct2), "data.frame")
    expect_equal(names(ct2), c("ME", "trait", "cor", "pvalue"))

})

