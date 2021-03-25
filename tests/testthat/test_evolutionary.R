
#----Load data----
data(og.zma.osa)
data(filt.se)

#----Start tests
test_that("is_singleton() returns logical vector of singletons", {
    genes <- tail(rownames(filt.se), n=100)
    s <- is_singleton(genes, og.zma.osa)
    expect_true(is.logical(s))
})
