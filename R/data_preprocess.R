#' Combine multiple expression tables (.tsv) into a single data frame
#'
#' This function reads multiple expression tables (.tsv files) in a directory and combine them into one single gene expression dataframe.
#'
#' @param mypath Path to directory containing .tsv files. Files must have the first column in common, e.g. "Gene_ID". Rows are gene IDs and columns are sample names.
#' @param pattern Pattern contained in each expression file. Default is '.tsv$', which means that all files ending in '.tsv' in the specified directory will be considered expression files.
#' @return Data frame with gene IDs as row names and their expression values in each sample (columns).
#' @author Fabricio Almeida-Silva
#' @rdname dfs2one
#' @export
#' @importFrom utils head read.csv write.table
dfs2one <- function(mypath, pattern = ".tsv$"){
    filenames <- list.files(path=mypath, full.names=TRUE, pattern = pattern)
    datalist <- lapply(filenames, function(x) {
        read.csv(file=x, header=T, sep="\t", stringsAsFactors=F, check.names=F)
    })
    merged.list <- Reduce(function(x,y) merge(x,y, all.x=T), datalist)
    rownames(merged.list) <- merged.list[,1]; merged.list[,1] <- NULL
    return(merged.list)
}

#' Remove missing values in a gene expression data frame
#'
#' @param exp A gene expression data frame with genes in row names and samples in column names or a `SummarizedExperiment` object.
#' @param replaceby What to use instead of NAs. One of 0 or 'mean'. Default is 0.
#'
#' @return Gene expression data frame or `SummarizedExperiment` object with all NAs replaced according to the argument 'replaceby'.
#' @author Fabricio Almeida-Silva
#' @export
#' @rdname remove_na
#' @examples
#' data(se.seed)
#' exp <- remove_na(se.seed)
#' sum(is.na(exp))
remove_na <- function(exp, replaceby = 0) {
    fexp <- handleSE(exp)

    if(is_SE(exp)) {
        fexp <- SummarizedExperiment::assay(exp)
    } else {
        fexp <- exp
    }

    if(replaceby == 0) {
        fexp[is.na(fexp)] <- 0
    } else {
        indices <- which(is.na(fexp), arr.ind = TRUE)
        fexp[indices] <- rowMeans(fexp, na.rm = TRUE)[indices[,1]]
    }

    if(is_SE(exp)) {
        fexp <- exp2SE(fexp, exp)
    }

    return(fexp)
}

#' Remove genes that are not expressed based on a user-defined threshold
#'
#' @param exp A gene expression data frame with genes in row names and samples in column names or a `SummarizedExperiment` object.
#' @param method Criterion to filter non-expressed genes out. One of "mean", "median", "percentage", or "allsamples". Default is "median".
#' @param min_exp If method is 'mean', 'median', or 'allsamples', the minimum value for a gene to be considered expressed. If method is 'percentage', the minimum value each gene must have in at least n percent of samples to be considered expressed.
#' @param min_percentage_samples In case the user chooses 'percentage' as method, expressed genes must have expression >= min_exp in at least this percentage. Values must range from 0 to 1.
#'
#' @return Filtered gene expression data frame or `SummarizedExperiment` object.
#' @author Fabricio Almeida-Silva
#' @export
#' @importFrom matrixStats rowMedians
#' @seealso
#'  \code{\link[matrixStats]{rowMedians}}
#'  \code{\link[WGCNA]{goodSamplesGenes}}
#' @rdname remove_nonexp
#' @examples
#' data(se.seed)
#' filt_exp <- remove_nonexp(se.seed, min_exp = 5)
remove_nonexp <- function(exp, method="median", min_exp=1, min_percentage_samples=0.25) {
    fexp <- handleSE(exp)

    if(method == "median") {
        final_exp <- fexp[matrixStats::rowMedians(as.matrix(fexp)) >= min_exp,]
    } else if (method == "mean") {
        final_exp <- fexp[rowMeans(fexp) >= min_exp,]
    } else if (method == "percentage") {
        min_n <- ncol(fexp) * min_percentage_samples
        final_exp <- fexp[rowSums(fexp >= min_exp) >= min_n, ]
    } else if (method == "allsamples") {
        final_exp <- fexp[rowSums(fexp >= min_exp) == ncol(fexp), ]
    } else {
        print("No method specified. Please, choose a filtering method - mean, median or percentage")
    }

    if(is_SE(exp)) {
        final_exp <- exp2SE(final_exp, exp)
    }

    return(final_exp)
}

#' Filter expression data frame to keep only the most varying genes
#'
#' @param exp A gene expression data frame with genes in row names and samples in column names or a `SummarizedExperiment` object.
#' @param n Number of most variable genes (e.g. n=5000 will keep the top 5000 most variable genes).
#' @param percentile Percentile of most highly variable genes (e.g. percentile=0.1 will keep the top 10 percent most variable genes). Values must range from 0 to 1.
#'
#' @return Gxpression data frame or `SummarizedExperiment` object with the most variable genes in row names and samples in column names.
#' @author Fabricio Almeida-Silva
#' @export
#' @rdname filter_by_variance
#' @examples
#' data(se.seed)
#' filt_exp <- filter_by_variance(se.seed, p=0.1)
filter_by_variance <- function(exp, n=NULL, percentile=NULL) {
    fexp <- handleSE(exp)

    gene_var <- data.frame(Genes = rownames(fexp), Var = apply(fexp, 1, var),
                           stringsAsFactors = FALSE)
    gene_var_ordered <- gene_var[order(gene_var$Var, decreasing = TRUE), ]
    if(!is.null(n) & is.null(percentile)) {
        top_var <- gene_var_ordered$Genes[seq_len(n)]
    } else if(is.null(n) & !is.null(percentile)) {
        p <- nrow(gene_var_ordered) * percentile
        top_var <- gene_var_ordered$Genes[seq_len(p)]
    } else {
        stop("Please, choose either 'n' or 'percentile'.")
    }
    top_variant_exp <- fexp[rownames(fexp) %in% top_var, ]

    if(is_SE(exp)) {
        top_variant_exp <- exp2SE(top_variant_exp, exp)
    }

    return(top_variant_exp)
}

#' Filter outlying samples based on standardized connectivity (Zk) method
#'
#' @param exp A gene expression data frame with genes in row names and samples in column names or a `SummarizedExperiment` object.
#' @param zk Standardized connectivity threshold. Default is -2.
#' @param cor_method Correlation method. One of "pearson", "biweight" or "spearman". Default is "spearman", considering that the expression data does not follow a normal distribution.
#'
#' @return Filtered gene expression data frame or `SummarizedExperiment` object.
#' @author Fabricio Almeida-Silva
#' @importFrom WGCNA adjacency bicor
#' @seealso
#'  \code{\link[WGCNA]{adjacency}}
#' @rdname ZKfiltering
#' @export
#' @examples
#' data(se.seed)
#' filt_exp <- ZKfiltering(se.seed)
#' @references
#' Oldham, M. C., Langfelder, P., & Horvath, S. (2012). Network methods for describing sample relationships in genomic datasets: application to Huntington’s disease. BMC systems biology, 6(1), 1-18.
ZKfiltering <- function(exp, zk = -2, cor_method = "spearman") {
    fexp <- handleSE(exp)

    if(cor_method == "pearson") {
        A <- WGCNA::adjacency(fexp, type = "distance")
    } else if(cor_method == "biweight") {
        A <- WGCNA::adjacency(fexp, type = "distance", corFnc = bicor, corOptions = list(use = 'p', maxPOutliers = 0.05))
    } else if(cor_method == "spearman"){
        A <- WGCNA::adjacency(fexp, type = "distance", corOptions = list(use = 'p', method = "spearman"))
    } else {
        print("Please, specify a correlation method (one of 'spearman', 'pearson' or 'biweight').")
    }

    k <- colSums(A) - 1
    Z.k <- scale(k)

    # Remove outliers
    remove.samples <- Z.k < zk | is.na(Z.k)
    fexp <- fexp[, !remove.samples]
    message("Number of removed samples: ", sum(remove.samples))

    if(is_SE(exp)) {
        fexp <- exp2SE(fexp, exp)
    }
    return(fexp)
}


#' Quantile normalize the expression data
#'
#' @param exp A gene expression data frame with genes in row names and samples in column names.
#'
#' @return Expression matrix with normalized values
#' @rdname q_normalize
#' @export
#' @examples
#' data(se.seed)
#' exp <- SummarizedExperiment::assay(se.seed)
#' norm_exp <- q_normalize(exp)
q_normalize <- function(exp) {
    n <- nrow(exp)
    rank.exp <- apply(exp, 2, rank)
    U <- rank.exp / (n+1)
    qnorm(U)
}

#' Apply Principal Component (PC)-based correction for confounding artifacts in coexpression networks
#'
#' @param exp A gene expression data frame with genes in row names and samples in column names or a `SummarizedExperiment` object.
#'
#' @return Corrected expression data frame or `SummarizedExperiment` object.
#' @author Fabricio Almeida-Silva
#' @export
#' @seealso
#'  \code{\link[sva]{num.sv}},\code{\link[sva]{sva_network}}
#' @rdname PC_correction
#' @importFrom sva num.sv sva_network
#' @examples
#' data(se.seed)
#' exp <- filter_by_variance(se.seed, n=500)
#' exp <- PC_correction(exp)
#' @references
#' Parsana, P., Ruberman, C., Jaffe, A. E., Schatz, M. C., Battle, A., & Leek, J. T. (2019). Addressing confounding artifacts in reconstruction of gene co-expression networks. Genome biology, 20(1), 1-6.
PC_correction <- function(exp) {
    fexp <- handleSE(exp)

    texp <- t(fexp) # transpose data frame
    texp <- as.matrix(texp)
    mod <- matrix(1, nrow = nrow(texp), ncol = 1)
    colnames(mod) <- "Intercept"

    message("Calculating number of PCs to be removed...")
    nsv <- sva::num.sv(t(texp), mod, method = "be")

    message("Number of PCs estimated to be removed: ", nsv)

    # PC residualization of gene expression data using sva_network
    message("Removing PCs that contribute to noise...")
    exprs_corrected <- sva::sva_network(texp, nsv)
    exprs_corrected_norm <- q_normalize(exprs_corrected)
    final.exp.corrected <- as.data.frame(exprs_corrected_norm)
    final.exp.corrected <- t(final.exp.corrected)

    if(is_SE(exp)) {
        final.exp.corrected <- exp2SE(final.exp.corrected, exp)
    }
    return(final.exp.corrected)
}


#' Preprocess expression data for network reconstruction
#'
#' @param exp A gene expression data frame with genes in row names and samples in column names or a `SummarizedExperiment` object.
#' @param NA_rm Logical. It specifies whether to remove missing values from the expression data frame or not. Default = TRUE.
#' @param replaceby If NA_rm is TRUE, what to use instead of NAs. One of 0 or 'mean'. Default is 0.
#' @param Zk_filtering Logical. It specifies whether to filter outlying samples by Zk or not. Default: TRUE.
#' @param zk If Zk_filtering is TRUE, the standardized connectivity threshold. Samples below this threshold will be considered outliers. Default is -2.
#' @param cor_method If Zk_filtering is TRUE, the correlation method to use. One of 'spearman', 'bicor', or 'pearson'. Default is 'spearman'.
#' @param remove_nonexpressed Logical. It specifies whether non-expressed genes should be removed or not. Default is TRUE.
#' @param method If remove_nonexpressed is TRUE, the criterion to filter non-expressed genes out. One of "mean", "median", "percentage", or "allsamples". Default is 'median'.
#' @param min_exp If method is 'mean', 'median', or 'allsamples', the minimum value for a gene to be considered expressed. If method is 'percentage', the minimum value each gene must have in at least n percent of samples to be considered expressed.
#' @param min_percentage_samples If method is 'percentage', expressed genes must have expression >= min_exp in at least this percentage. Values must range from 0 to 1. Default = 0.25.
#' @param remove_confounders Logical. If TRUE, it removes principal components that add noise to the data.
#' @param variance_filter Logical. If TRUE, it will filter genes by variance. Default is FALSE.
#' @param n If variance_filter is TRUE, the number of the most variable genes to keep.
#' @param percentile If variance_filter is TRUE, the percentage of the most variable genes to keep.
#' @param vstransform Logical indicating if data should be variance stabilizing transformed. This parameter can only be set to TRUE if data is a matrix of raw read counts.
#'
#' @return Processed gene expression data frame or `SummarizedExperiment` object with gene IDs in row names and sample names in column names.
#' @author Fabricio Almeida-Silva
#' @seealso
#'  \code{\link[DESeq2]{varianceStabilizingTransformation}}
#' @rdname exp_preprocess
#' @export
#' @importFrom DESeq2 varianceStabilizingTransformation
#' @examples
#' data(se.seed)
#' exp <- exp_preprocess(se.seed, variance_filter=TRUE, n=1000)
#' @references
#' Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome biology, 15(12), 1-21.
exp_preprocess <- function(exp, NA_rm = TRUE, replaceby = 0,
                           Zk_filtering = TRUE, zk = -2,
                           cor_method = "spearman",
                           remove_nonexpressed = TRUE, method = "median",
                           min_exp = 1, min_percentage_samples = 0.25,
                           remove_confounders = TRUE,
                           variance_filter = FALSE, n = NULL, percentile = NULL,
                           vstransform = FALSE) {
    # Remove missing values
    if(NA_rm) {
        exp <- remove_na(exp, replaceby = replaceby)
    }
    # Remove non-expressed genes
    if(remove_nonexpressed) {
        exp <- remove_nonexp(exp, method = method,
                             min_exp = min_exp,
                             min_percentage_samples = min_percentage_samples)
    }
    # Apply VST for count data
    if(vstransform) {
        fexp <- handleSE(exp)
        fexp <- as.data.frame(DESeq2::varianceStabilizingTransformation(fexp))
        if(is_SE(exp)) {
            fexp <- exp2SE(fexp, exp)
        }
        exp <- fexp
    }
    # Filter by variance
    if(variance_filter) {
        exp <- filter_by_variance(exp, n = n, percentile = percentile)
    }
    # Zk filtering
    if(Zk_filtering) {
        exp <- ZKfiltering(exp, zk = zk, cor_method = cor_method)
    }
    # Remove confounders
    if(remove_confounders) {
        exp <- PC_correction(exp)
    }
    return(exp)
}


#' Create an RMA-normalized expression data frame from .cel files
#'
#' @param cel_dir Path to directory containing binary .cel files. Default is the current working directory.
#' @param pattern Character indicating pattern to match binary .cel files. Default is "cel", indicating that all .cel files must have "cel" in their file names.
#'
#' @return An RMA-normalized expression data frame with probe IDs as row names and arrays in column names.
#' @seealso
#'  \code{\link[affy]{read.affybatch}},\code{\link[affy]{rma}}
#'  \code{\link[Biobase]{eSet}},\code{\link[Biobase]{exprs}}
#' @rdname cel2RMA
#' @export
#' @importFrom affy ReadAffy rma
#' @importFrom Biobase exprs
cel2RMA <- function(cel_dir = "./", pattern = ".cel") {
    cels <- list.files(cel_dir, pattern = pattern)
    raw.data <- affy::ReadAffy(verbose = FALSE, filenames = cels)
    rma.norm <- affy::rma(raw.data)
    rma.exp <- Biobase::exprs(rma.norm)

    colnames(rma.exp) <- sub(".cel", "", colnames(rma.exp), ignore.case = TRUE)

    return(rma.exp)
}

#' Collapse probe-level expression matrix (or data frame) to gene-level
#'
#' @param probe_exp Expression matrix or data frame containing probe names as row names and chips as column names
#' @param correspondence A 2-column data frame containing gene IDs in the first column and their corresponding probe names in the second column. Probe names must be unique, but gene IDs can be repeated, since different probes can match the same gene.
#' @return A data frame containing gene IDs as row names and arrays in the column names.
#' @details If there are 2 or more probes for the same gene, the one with the highest mean will be selected.
#' @seealso
#'  \code{\link[WGCNA]{collapseRows}}
#' @rdname probe2gene
#' @export
#' @importFrom WGCNA collapseRows
probe2gene <- function(probe_exp, correspondence) {
    colnames(correspondence) <- c("rowGroup", "rowID")
    correspondence <- correspondence[correspondence$rowID %in% rownames(probe_exp), ]

    exp <- WGCNA::collapseRows(probe_exp, rowGroup = correspondence$rowGroup,
                               rowID = correspondence$rowID)

    final_exp <- as.data.frame(exp[[1]])
    return(final_exp)
}
