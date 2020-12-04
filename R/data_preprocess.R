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
#' @param exp Gene expression data frame with genes in row names and samples IDs in column names.
#' @param replaceby What to use instead of NAs. One of 0 or 'mean'. Default is 0.
#'
#' @return Gene expression data frame with all NAs replaced according to the argument 'replaceby'
#' @author Fabricio Almeida-Silva
#' @export
#' @rdname remove_na
#' @export

remove_na <- function(exp, replaceby = 0) {
    if(replaceby == 0) {
        exp[is.na(exp)] <- 0
    } else {
        indices <- which(is.na(exp), arr.ind = TRUE)
        exp[indices] <- rowMeans(exp, na.rm = TRUE)[indices[,1]]
    }
    return(exp)
}

#' Remove genes that are not expressed based on a user-defined threshold
#'
#' @param exp Gene expression data frame, where rownames are gene IDs and colnames are sample names.
#' @param method Criterion to filter non-expressed genes out. One of "mean", "median", "percentage", or "allsamples". Default is "median".
#' @param min_exp If method is 'mean', 'median', or 'allsamples', the minimum value for a gene to be considered expressed. If method is 'percentage', the minimum value each gene must have in at least n percent of samples to be considered expressed.
#' @param min_percentage_samples in case the user chooses 'percentage' as method, expressed genes must have expression >= min_exp in at least this percentage. Values must range from 0 to 1.
#'
#' @return Dataframe containing only expressed genes in rownames and sample names in colnames
#' @author Fabricio Almeida-Silva
#' @export
#' @importFrom matrixStats rowMedians
#' @importFrom WGCNA goodSamplesGenes
#' @importFrom dynamicTreeCut printFlush
#' @seealso
#'  \code{\link[matrixStats]{rowMedians}}
#'  \code{\link[WGCNA]{goodSamplesGenes}}
#' @rdname remove_nonexp
remove_nonexp <- function(exp, method="median", min_exp=1, min_percentage_samples=0.25) {
    if(method == "median") {
        final_exp <- exp[matrixStats::rowMedians(as.matrix(exp)) >= min_exp,]
    } else if (method == "mean") {
        final_exp <- exp[rowMeans(exp) >= min_exp,]
    } else if (method == "percentage") {
        exp[exp < min_exp] <- NA
        texp <- t(exp)

        gsg = WGCNA::goodSamplesGenes(texp, minFraction = min_percentage_samples, verbose = 3)
        if (!gsg$allOK)
        {
            if (sum(!gsg$goodGenes)>0)
                dynamicTreeCut::printFlush(paste("Removing genes:", paste(names(texp)[!gsg$goodGenes], collapse = "\t")))
            if (sum(!gsg$goodSamples)>0)
                dynamicTreeCut::printFlush(paste("Removing samples:", paste(rownames(texp)[!gsg$goodSamples], collapse = "\t")))
            final_exp = exp[gsg$goodGenes, gsg$goodSamples]
        }
    } else if (method == "allsamples") {
        final_exp <- exp[rowSums(exp >= min_exp) == ncol(exp), ]
    } else {
        print("No method specified. Please, choose a filtering method - mean, median or percentage")
    }
    final_exp[is.na(final_exp)] <- 0
    return(final_exp)
}

#' Filter expression data frame to keep only the most varying genes
#'
#' @param exp Expression data frame with genes as rows and samples as columns
#' @param n Number of most variable genes (e.g. n=5000 will keep the top 5000 most variable genes).
#' @param percentile Percentile of most highly variable genes (e.g. percentile=0.1 will keep the top 10 percent most variable genes). Values must range from 0 to 1.
#'
#' @return Expressed data frame with the most variable genes in row names and samples in column names.
#' @author Fabricio Almeida-Silva
#' @export
#' @rdname filter_by_variance

filter_by_variance <- function(exp, n=NULL, percentile=NULL) {
    gene_var <- data.frame(Genes = rownames(exp), Var = apply(exp, 1, var), stringsAsFactors = FALSE)
    gene_var_ordered <- gene_var[order(gene_var$Var, decreasing = TRUE), ]
    if(!is.null(n) & is.null(percentile)) {
        top_var <- gene_var_ordered$Genes[1:n]
    } else if(is.null(n) & !is.null(percentile)) {
        p <- nrow(gene_var_ordered) * percentile
        top_var <- gene_var_ordered$Genes[1:p]
    } else {
        stop("Please, choose either 'n' or 'percentile'.")
    }
    top_variant_exp <- exp[rownames(exp) %in% top_var, ]
}

#' Filter outlying samples based on standardized connectivity (Zk) method
#'
#' @param raw_exp Raw gene expression table, where rownames are gene IDs and colnames are sample names.
#' @param Zk Standardized connectivity threshold. Default is -2.
#' @param cor_method Correlation method. One of "pearson", "biweight" or "spearman". Default is "spearman", considering that the expression data does not follow a normal distribution.
#'
#' @return Filtered gene expression dataframe.
#' @author Fabricio Almeida-Silva
#' @importFrom WGCNA adjacency bicor
#' @seealso
#'  \code{\link[WGCNA]{adjacency}}
#' @rdname ZKfiltering
#' @export
ZKfiltering <- function(raw_exp, Zk = -2, cor_method = "spearman") {
    if(cor_method == "pearson") {
        A = WGCNA::adjacency(raw_exp, type = "distance")
    } else if(cor_method == "biweight") {
        A = WGCNA::adjacency(raw_exp, type = "distance", corFnc = bicor, corOptions = list(use = 'p', maxPOutliers = 0.05))
    } else if(cor_method == "spearman"){
        A = WGCNA::adjacency(raw_exp, type = "distance", corOptions = list(use = 'p', method = "spearman"))
    } else {
        print("Please, specify a correlation method (one of 'spearman', 'pearson' or 'biweight').")
    }

    k = as.numeric(apply(A, 2, sum))-1
    Z.k = scale(k)

    #Say that every sample whose Z.k value is below <Zk> is an outlier
    thresholdZ.k = Zk

    #Remove outliers
    remove.samples <- Z.k < thresholdZ.k | is.na(Z.k)
    traw_exp <- t(raw_exp)
    traw_exp2 <- traw_exp[!remove.samples,]

    #Create a new dataframe with filtered samples and genes
    print(paste("Number of samples that were removed:", sum(remove.samples)))
    raw_exp_filtsamp <- t(as.data.frame(traw_exp2))
    return(raw_exp_filtsamp)
}

#' Preprocess expression data for network reconstruction
#'
#' @param exp Gene Expression data frame with gene IDs as rownames and sample names as column names.
#' @param NA_rm Logical. It specifies whether to remove missing values from the expression data frame or not. Default = TRUE.
#' @param replaceby If NA_rm is TRUE, what to use instead of NAs. One of 0 or 'mean'. Default is 0.
#' @param Zk_filtering Logical. It specifies whether to filter outlying samples by Zk or not. Default = TRUE.
#' @param Zk If Zk_filtering is TRUE, the standardized connectivity threshold. Samples below this threshold will be considered outliers. Default is -2.
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
#' @return Processed gene expression data frame with gene IDs in row names and sample names in column names.
#' @author Fabricio Almeida-Silva
#' @seealso
#'  \code{\link[DESeq2]{varianceStabilizingTransformation}}
#' @rdname exp_preprocess
#' @export
#' @importFrom DESeq2 varianceStabilizingTransformation
exp_preprocess <- function(exp, NA_rm = TRUE, replaceby = 0, Zk_filtering = TRUE, Zk = -2, cor_method = "spearman",
                           remove_nonexpressed = TRUE, method = "median", min_exp = 1, min_percentage_samples = 0.25,
                           remove_confounders = TRUE, variance_filter = FALSE, n = NULL, percentile = NULL,
                           vstransform = FALSE) {

    # Remove missing values
    if(NA_rm == TRUE) {
        exp1 <- remove_na(exp, replaceby = replaceby)
    } else {
        exp1 <- exp
    }

    # Remove non-expressed genes
    if(remove_nonexpressed == TRUE) {
        exp2 <- remove_nonexp(exp1, method = method, min_exp = min_exp, min_percentage_samples = min_percentage_samples)
    } else {
        exp2 <- exp1
    }

    if(vstransform == TRUE) {
        exp2 <- as.data.frame(DESeq2::varianceStabilizingTransformation(exp2))
    }

    # Filter by variance
    if(variance_filter == TRUE) {
        exp3 <- filter_by_variance(exp2, n = n, percentile = percentile)
    } else {
        exp3 <- exp2
    }

    # Zk filtering
    if(Zk_filtering == TRUE) {
        exp4 <- ZKfiltering(exp3, Zk = Zk, cor_method = cor_method)
    } else {
        exp4 <- exp3
    }

    # Remove confounders
    if(remove_confounders == TRUE) {
        exp5 <- as.data.frame(PC_correction(exp4))
    } else {
        exp5 <- exp4
    }

    return(exp5)
}

#' Quantile normalize the expression data
#'
#' @param exp Gene expression matrix, where rownames are sample names and colnames are gene IDs.
#'
#' @return Expression matrix with normalized values
#' @rdname q_normalize
#' @export

q_normalize <- function(exp){
    n <- nrow(exp)
    p <- ncol(exp)
    rank.exp <- exp # matrix for ranking
    for (i in 1:p){
        rank.exp[,i] <- rank(exp[,i])
    }
    U <- rank.exp/(n+1)
    qnorm(U)
}

#' Apply Principal Component (PC)-based correction for confounding artifacts in coexpression networks
#'
#' @param exp_df Expression dataframe where rownames are gene IDs and colnames are sample names. Expression data must not have undergone normalization.
#'
#' @return Corrected expression data frame
#' @author Fabricio Almeida-Silva
#' @export
#' @seealso
#'  \code{\link[sva]{num.sv}},\code{\link[sva]{sva_network}}
#' @rdname PC_correction
#' @importFrom sva num.sv sva_network
PC_correction <- function(exp_df) {
    texp <- t(exp_df) #transpose data so that rows are samples and columns are genes
    texp <- as.matrix(texp)
    mod <- matrix(1, nrow=dim(texp)[1], ncol=1)
    colnames(mod) <- "Intercept"

    print("Calculating number of PCs to be removed...")
    nsv <- sva::num.sv(t(texp), mod, method = "be") # num.sv requires data matrix with features(genes) in the rows and samples in the column

    print(paste("Number of PCs estimated to be removed:", nsv))

    #PC residualization of gene expression data using sva_network
    print("Removing PCs that contribute to noise...")
    exprs_corrected <- sva::sva_network(texp, nsv)
    exprs_corrected_norm <- q_normalize(exprs_corrected)
    final.exp.corrected <- as.data.frame(exprs_corrected_norm)
    final.exp.corrected <- t(final.exp.corrected)
    return(final.exp.corrected)
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
cel2RMA <- function(cel_dir = "./", pattern = "cel") {
    cels <- list.files(cel_dir, pattern = pattern)
    raw.data <- affy::ReadAffy(verbose = FALSE, filenames = cels)
    rma.norm <- affy::rma(raw.data)
    rma.exp <- Biobase::exprs(rma.norm)

    colnames(rma.exp) <- sub(".cel", "", colnames(rma.exp), ignore.case = TRUE)

    return(rma.exp)
}

#' Collapse probe-level expression matrix (or data frame) to gene-level
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
