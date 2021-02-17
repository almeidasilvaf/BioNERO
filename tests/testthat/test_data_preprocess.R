
#----Set up data sets----
data(zma.se)
exp <- SummarizedExperiment::assay(zma.se)
metadata <- SummarizedExperiment::colData(zma.se)
metadata <- data.frame(Samples=rownames(metadata),
                       Tissue=metadata$Tissue)

#----Start tests----
test_that("dfs2one reads in multiple tables and binds them into a data frame", {
    genes <- paste0(rep("Gene", 100), 1:100)
    samples1 <- paste0(rep("Sample", 30), 1:30)
    samples2 <- paste0(rep("Sample", 30), 31:60)

    # Simulate two expression data frames of 100 genes and 30 samples
    exp1 <- cbind(genes, as.data.frame(matrix(rnorm(100*30),nrow=100,ncol=30)))
    exp2 <- cbind(genes, as.data.frame(matrix(rnorm(100*30),nrow=100,ncol=30)))
    colnames(exp1) <- c("Gene", samples1)
    colnames(exp2) <- c("Gene", samples2)

    # Write data frames to temporary files
    tmpdir <- tempdir()
    tmp1 <- tempfile(tmpdir = tmpdir, fileext = ".exp.tsv")
    tmp2 <- tempfile(tmpdir = tmpdir, fileext = ".exp.tsv")

    write.table(exp1, file=tmp1, quote=FALSE, sep="\t")
    write.table(exp2, file=tmp2, quote=FALSE, sep="\t")

    expect_equal(nrow(dfs2one(mypath = tmpdir, pattern=".exp.tsv")), 100)
    expect_equal(ncol(dfs2one(mypath = tmpdir, pattern=".exp.tsv")), 60)
})


test_that("remove_na() removes NAs for both SE and expression data frame", {
    expinput <- remove_na(exp)
    seinput <- remove_na(zma.se)
    expect_equal(sum(is.na(expinput)), 0)
    expect_equal(sum(is.na(seinput)), 0)
    expect_equal(dim(expinput), dim(exp))
    expect_equal(dim(seinput), dim(zma.se))
})


test_that("remove_nonexp() removes non-expressed genes", {
    filt_exp1 <- remove_nonexp(exp, min_exp=10)
    filt_exp2 <- remove_nonexp(zma.se, min_exp=10)
    filt_exp3 <- remove_nonexp(zma.se, method = "percentage", min_exp=10)
    filt_exp4 <- remove_nonexp(zma.se, method="mean", min_exp=10)
    filt_exp5 <- remove_nonexp(zma.se, method="allsamples", min_exp=10)

    expect_true(all.equal(dim(filt_exp1), dim(filt_exp2)))
    expect_equal(class(filt_exp1), "data.frame")
    expect_true(class(filt_exp2) == "SummarizedExperiment")
    expect_true(nrow(filt_exp2) < nrow(zma.se))
    expect_true(nrow(filt_exp3) < nrow(zma.se))
    expect_true(nrow(filt_exp4) < nrow(zma.se))
    expect_true(nrow(filt_exp5) < nrow(zma.se))
})


test_that("filter_by_variance() filters gene expression data by variance", {
    filt_exp1 <- filter_by_variance(exp, n=3000)
    filt_exp2 <- filter_by_variance(zma.se, n=3000)

    expect_true(all.equal(dim(filt_exp1), dim(filt_exp2)))
    expect_true(class(filt_exp2) == "SummarizedExperiment")
    expect_equal(class(filt_exp1), "data.frame")
    expect_equal(nrow(filt_exp1), 3000)
    expect_equal(nrow(filt_exp2), 3000)
})


test_that("ZKfiltering() filters outliers", {
    filt_exp1 <- ZKfiltering(exp)
    filt_exp2 <- ZKfiltering(zma.se)
    filt_exp3 <- ZKfiltering(zma.se, cor_method = "pearson")
    filt_exp4 <- ZKfiltering(zma.se, cor_method = "biweight", zk = -2.5)

    expect_true(class(filt_exp2) == "SummarizedExperiment")
    expect_true(all.equal(dim(filt_exp1), dim(filt_exp2)))
    expect_message(filt_exp3 <- ZKfiltering(zma.se, cor_method = "pearson"),
                   "Number of removed samples: 1"
    )
})


test_that("q_normalize() performs quantile normalization in an expression data frame", {
    norm_exp <- q_normalize(exp)
    cols <- sample(seq_len(ncol(norm_exp)), size=2, replace=FALSE)

    expect_equal(
        IQR(norm_exp[, cols[1]]), IQR(norm_exp[, cols[2]])
    )
})


test_that("PC_correction() removes confounders based on principal components", {
    filt_exp <- filter_by_variance(exp, n=500)
    filt_exp2 <- filter_by_variance(zma.se, n=500)
    pc_cor1 <- PC_correction(filt_exp)
    pc_cor2 <- PC_correction(filt_exp2)

    expect_equal(dim(pc_cor1), dim(pc_cor2))
})


test_that("exp_preprocess() preprocesses the expression data at once", {
    filt_exp1 <- exp_preprocess(exp, variance_filter = TRUE, n=500)
    filt_exp2 <- exp_preprocess(zma.se, variance_filter = TRUE, n=500)

    expect_equal(dim(filt_exp1), dim(filt_exp2))
})









