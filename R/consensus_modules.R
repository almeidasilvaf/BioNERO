

#' Pick power to fit networks to scale-free topology
#'
#' @param exp_list A list of expression data frames or
#' SummarizedExperiment objects.
#' If input is a list of data frames, row names must correspond to gene IDs
#' and column names to samples.
#' The list can be created with \code{list(exp1, exp2, ..., expn)}.
#' @param setLabels Character vector containing labels for each expression set.
#' @param metadata A data frame containing sample names in row names and
#' sample annotation in the first column.
#' Ignored if `exp_list` is a list of `SummarizedExperiment` objects, since
#' the function will extract colData.
#' @param cor_method Correlation method used for network reconstruction.
#' One of "spearman" (default), "biweight", or "pearson".
#' @param net_type Network type. One of "signed hybrid" (default),
#' "signed" or "unsigned".
#' @param rsquared Minimum R squared to consider the network similar to
#' a scale-free topology. Default is 0.8.
#'
#' @return A list of 2 elements: \describe{
#'   \item{power}{Numeric vector of optimal beta powers to fit networks to SFT}
#'   \item{plot}{A ggplot object displaying main statistics of the SFT fit test}
#' }
#' @rdname consensus_SFT_fit
#' @export
#' @importFrom WGCNA pickSoftThreshold checkSets
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw scale_color_manual
#' ylim theme geom_hline
#' @importFrom ggrepel geom_text_repel
#' @examples
#' set.seed(12)
#' data(zma.se)
#' filt.zma <- filter_by_variance(zma.se, n=500)
#' zma.set1 <- filt.zma[, sample(colnames(filt.zma), size=20, replace=FALSE)]
#' zma.set2 <- filt.zma[, sample(colnames(filt.zma), size=20, replace=FALSE)]
#' list.sets <- list(zma.set1, zma.set2)
#' cons_sft <- consensus_SFT_fit(list.sets, setLabels = c("Maize1", "Maize2"),
#'                               cor_method = "pearson")
consensus_SFT_fit <- function(exp_list, setLabels = NULL, metadata = NULL,
                              cor_method = "spearman",
                              net_type = "signed hybrid", rsquared = 0.8) {
    metadata <- handle_metadata(exp_list, metadata)
    exp_list <- handleSElist(exp_list)
    if(is.null(setLabels)) {
        setLabels <- seq_along(exp_list)
    }

    # Build multi-set object
    multiExp <- lapply(exp_list, function(x) {
        element <- list(data=as.data.frame(t(x)))
        return(element)
    })

    # Check if the data has the correct format for downstream analysis
    expSize <- WGCNA::checkSets(multiExp)

    # Choose a set of soft-thresholding powers
    powers <- seq(5, 20, by=1)
    sft <- lapply(multiExp, function(x) {
        if(cor_method == "pearson") {
            sft <- WGCNA::pickSoftThreshold(
                x$data, networkType = net_type, powerVector = powers,
                RsquaredCut = rsquared
            )
        } else if(cor_method == "biweight") {
            sft <- WGCNA::pickSoftThreshold(
                x$data, networkType = net_type, powerVector = powers,
                RsquaredCut = rsquared, corFnc = bicor,
                corOptions = list(use = 'p', maxPOutliers = 0.05)
            )
        } else if (cor_method == "spearman") {
            sft <- WGCNA::pickSoftThreshold(
                x$data, networkType = net_type, powerVector = powers,
                RsquaredCut = rsquared,
                corOptions = list(use = 'p', method = "spearman")
            )
        } else {
            stop("Please, specify a correlation method (one of 'spearman', 'pearson' or 'biweight').")
        }
        wgcna_power <- sft$powerEstimate
        if(is.na(wgcna_power)){
            message("No power reached R-squared cut-off, now choosing max R-squared based power")
            wgcna_power <- sft$fitIndices$Power[which(sft$fitIndices$SFT.R.sq == max(sft$fitIndices$SFT.R.sq))]
        }
        return(list(sft = sft, wgcna_power = wgcna_power))
    })

    power <- unlist(lapply(sft, function(x) return(x$wgcna_power)))

    # Create data frame with indices to plot
    sft_df <- Reduce(rbind, lapply(seq_along(sft), function(x) {
        indices <- sft[[x]]$sft$fitIndices
        df <- data.frame(power = indices[,1],
                         fit = -sign(indices[,3]) * indices[,2],
                         meank = indices[,5],
                         Set = setLabels[x])
        return(df)
    }))

    # Plot 1
    cols <- custom_palette(1)
    p1 <- ggplot(sft_df, aes(x = .data$power, y = .data$fit)) +
        geom_point(aes(color = .data$Set)) +
        ggrepel::geom_text_repel(aes(label = .data$power, color = .data$Set)) +
        scale_color_manual(values = cols) +
        labs(
            x = "Soft threshold (power)",
            y = expression(paste("Scale-free topology fit - ", R^{2})),
            title = "Scale independence"
        ) +
        theme_bw() +
        ylim(c(0,1)) +
        geom_hline(yintercept = rsquared, color = "brown3") +
        theme(legend.position = "none")

    # Plot 2
    p2 <- ggplot(sft_df, aes(x = .data$power, y = .data$meank)) +
        geom_point(aes(color = .data$Set)) +
        ggrepel::geom_text_repel(aes(label = .data$power, color = .data$Set),
                                 show.legend = FALSE) +
        scale_color_manual(values = cols) +
        labs(
            x = "Soft threshold (power)",
            y = "Mean connectivity (k)",
            title = "Mean connectivity"
        ) +
        theme_bw()

    sft_plot <- patchwork::wrap_plots(p1, p2, ncol = 2)

    results <- list(power = power, plot = sft_plot)
    return(results)
}


#' Identify consensus modules across independent data sets
#'
#' @param exp_list A list containing the expression data frames with genes in
#' row names and samples in column names or `SummarizedExperiment` objects.
#' The list can be created by using \code{list(exp1, exp2, ..., expn)}.
#' @param metadata A data frame containing sample names in row names and
#' sample annotation in the first column.
#' Ignored if `exp_list` is a list of `SummarizedExperiment` objects, since
#' the function will extract colData.
#' @param power Numeric vector of beta power for each expression set
#' as calculated by \code{consensus_SFT_fit}.
#' @param cor_method Correlation method used for network reconstruction.
#' One of "spearman" (default), "biweight", or "pearson".
#' @param net_type Network type. One of "signed hybrid" (default),
#' "signed" or "unsigned".
#' @param module_merging_threshold Correlation threshold to merge
#' similar modules into a single one. Default: 0.8.
#' @param TOM_type Character indicating the type of Topological
#' Overlap Matrix to (TOM) create. One of
#' 'unsigned', 'signed', 'signed Nowick', 'unsigned 2', 'signed 2',
#' and 'signed Nowick 2'. By default, TOM type is automatically
#' selected based on network type.
#' @param verbose Logical indicating whether to display progress
#' messages or not. Default: FALSE.
#'
#' @return A list containing 4 elements: \describe{
#'   \item{consMEs}{Consensus module eigengenes}
#'   \item{exprSize}{Description of the multi-set object returned by the function \code{WGCNA::checkSets}}
#'   \item{sampleInfo}{Metadata for each expression set}
#'   \item{genes_cmodules}{Data frame of genes and consensus modules}
#'   \item{dendro_plot_objects}{Objects to be used in dendrogram plotting}
#' }
#'
#' @rdname consensus_modules
#' @export
#' @importFrom WGCNA checkSets adjacency TOMsimilarity pquantile labels2colors
#' multiSetMEs consensusMEDissimilarity mergeCloseModules
#' @importFrom dynamicTreeCut cutreeDynamic
#' @examples
#' set.seed(12)
#' data(zma.se)
#' filt.zma <- filter_by_variance(zma.se, n=500)
#' zma.set1 <- filt.zma[, sample(colnames(filt.zma), size=20, replace=FALSE)]
#' zma.set2 <- filt.zma[, sample(colnames(filt.zma), size=20, replace=FALSE)]
#' list.sets <- list(zma.set1, zma.set2)
#' # SFT power previously identified with consensus_SFT_fit()
#' cons_mod <- consensus_modules(list.sets, power = c(11, 13),
#'                               cor_method = "pearson")
consensus_modules <- function(
        exp_list, metadata, power,
        cor_method = "spearman", net_type = "signed hybrid",
        module_merging_threshold = 0.8,
        TOM_type = NULL,
        verbose = FALSE
) {

    nSets <- length(exp_list)

    # Extract sample metadata as a list of data frames
    if(is(exp_list[[1]], "SummarizedExperiment")) {
        metadata <- lapply(exp_list, function(x) {
            return(se2metadata(x)$coldata)
        })
    }

    # Keep only shared columns in sample metadata data frames
    shared <- Reduce(intersect, lapply(metadata, colnames))
    metadata <- lapply(metadata, function(x) return(x[, shared, drop = FALSE]))

    # Build multi-set object
    exp_list <- handleSElist(exp_list)
    multiExp <- lapply(exp_list, function(x) {
        element <- list(data = as.data.frame(t(x)))
        return(element)
    })

    # Check if the data has the correct format for downstream analysis
    expSize <- WGCNA::checkSets(multiExp)
    nGenes <- expSize$nGenes
    nSamples <- expSize$nSamples

    # Build multi-set object with metadata
    sampleinfo <- lapply(seq_along(multiExp), function(x) {
        sinfo <- metadata[[x]][rownames(multiExp[[x]]$data), , drop = FALSE]
        return(sinfo)
    })

    # Calculate adjacencies for each individual set
    if(verbose) { message("Calculating adjacency matrix...") }
    adj <- lapply(seq_len(nSets), function(x) {
        if(cor_method == "pearson") {
            adjacencies <- WGCNA::adjacency(
                multiExp[[x]]$data, power = power[x], type = net_type
            )
        } else if(cor_method == "spearman") {
            adjacencies <- WGCNA::adjacency(
                multiExp[[x]]$data, power = power[x], type = net_type,
                corOptions = list(use = "p", method = "spearman")
            )
        } else if (cor_method == "biweight") {
            adjacencies <- WGCNA::adjacency(
                multiExp[[x]]$data, power = power[x], type = net_type,
                corFnc = bicor
            )
        } else {
            stop("Please, specify a correlation method. One of 'spearman', 'pearson' or 'biweight'.")
        }
        return(adjacencies)
    })

    # Calculate TOMs in each individual data set
    if(verbose) { message("Calculating topological overlap matrix (TOM)...") }
    tom_type <- switch(
        net_type,
        "signed hybrid" = "signed",
        "signed" = "signed Nowick",
        "unsigned"
    )
    if(!is.null(TOM_type)) { tom_type <- TOM_type }

    TOM <- lapply(adj, function(x) return(TOMsimilarity(x, TOMType = tom_type)))

    # Scaling TOM to make them comparable
    scaleP <- 0.95
    nSamples <- as.integer(1 / (1-scaleP) * 100)
    scaleSample <- sample(nGenes * (nGenes-1) / 2, size = nSamples)

    scaledTOM <- lapply(seq_len(nSets), function(x) {
        tom_scaled <- TOM[[x]]
        if(x > 1) {
            TOMScalingSamples <- as.dist(TOM[[x]])[scaleSample]
            scaleQuant <- quantile(TOMScalingSamples, probs = scaleP, type = 8)
            scalePowers <- log(scaleQuant[1]) / log(scaleQuant)
            tom_scaled <- tom_scaled^scalePowers
        }
        return(tom_scaled)
    })

    # Calculation of consensus TOM
    if (nSets <= 3) {
        consensusTOM <- do.call(pmin, lapply(seq_along(TOM), function(i) TOM[[i]]))
    } else {
        consensusTOM <- do.call(function(x) WGCNA::pquantile(x, prob=0.25),
                                lapply(seq_along(TOM), function(i) TOM[[i]]))
    }

    # Clustering and module identification
    consTree <- hclust(as.dist(1-consensusTOM), method = "average")
    unmergedLabels <- dynamicTreeCut::cutreeDynamic(
        dendro = consTree,
        distM = 1 - consensusTOM,
        deepSplit = 2,
        cutHeight = 0.995,
        minClusterSize = 30,
        pamRespectsDendro = FALSE
    )

    nmod <- length(unique(unmergedLabels))
    palette <- rev(WGCNA::standardColors(nmod))
    unmergedColors <- WGCNA::labels2colors(unmergedLabels, colorSeq = palette)

    # Merge similar modules
    unmergedMEs <- WGCNA::multiSetMEs(
        multiExp, colors = NULL, universalColors = unmergedColors
    )

    merging_threshold <- 1-module_merging_threshold
    merge <- WGCNA::mergeCloseModules(
        multiExp, unmergedColors, cutHeight = merging_threshold, verbose = 0
    )

    genes_cmod <- data.frame(
        Genes = rownames(adj[[1]]),
        Cons_modules = merge$colors
    )
    consMEs <- merge$newMEs

    # Plot dendrogram with merged colors
    result_list <- list(
        consMEs = consMEs,
        exprSize = expSize,
        sampleInfo = sampleinfo,
        genes_cmodules = genes_cmod,
        dendro_plot_objects = list(
            tree = consTree,
            Unmerged = unmergedColors,
            Merged = moduleColors
        )
    )
    return(result_list)
}


#' Correlate set-specific modules and consensus modules to sample information
#'
#' @param consensus Consensus network returned by \code{consensus_modules}.
#' @param cor_method Correlation method to be used. One of 'spearman' or
#' 'pearson'. Default: 'pearson'.
#' @param metadata_cols A vector (either numeric or character) indicating
#' which columns should be extracted from column metadata if \strong{exp}
#' is a `SummarizedExperiment` object. The vector can contain column
#' indices (numeric) or column names (character). By default, all columns are
#' used.
#'
#' @return Data frame of consensus module-trait correlations and p-values,
#' with the following variables:
#' \describe{
#'   \item{trait}{Factor, trait name. Each trait corresponds to a variable
#'                of the sample metadata (if numeric) or levels of a variable
#'                (if categorical).}
#'   \item{ME}{Factor, module eigengene.}
#'   \item{cor}{Numeric, correlation.}
#'   \item{pvalue}{Numeric, correlation P-values.}
#'   \item{group}{Character, name of the metadata variable.}
#' }
#'
#' @rdname consensus_trait_cor
#' @export
#' @importFrom WGCNA corPvalueFisher labels2colors
#' @importFrom reshape2 melt
#' @examples
#' set.seed(12)
#' data(zma.se)
#' filt.zma <- filter_by_variance(zma.se, n=500)
#' zma.set1 <- filt.zma[, sample(colnames(filt.zma), size=20, replace=FALSE)]
#' zma.set2 <- filt.zma[, sample(colnames(filt.zma), size=20, replace=FALSE)]
#' list.sets <- list(zma.set1, zma.set2)
#' # SFT power previously identified with consensus_SFT_fit()
#' consensus <- consensus_modules(list.sets, power = c(11, 13),
#'                                cor_method = "pearson")
#' consensus_trait <- consensus_trait_cor(consensus, cor_method = "pearson")
consensus_trait_cor <- function(consensus, cor_method = "pearson",
                                metadata_cols = NULL) {

    consMEs <- consensus$consMEs
    metadata <- consensus$sampleInfo
    nSets <- consensus$exprSize$nSets

    if(!is.null(metadata_cols)) {
        metadata <- lapply(metadata, function(x) return(x[, metadata_cols, drop = FALSE]))
    }


    # Get a list of correlation and P-value matrices
    matrices <- lapply(seq_len(nSets), function(set) {

        coldata <- metadata[[set]]

        # Get model matrices
        mats <- Reduce(cbind, lapply(seq_along(coldata), function(x) {
            model_mat <- get_model_matrix(coldata, column_idx = x)
            return(model_mat)
        }))

        # Get correlation matrix
        cormat <- cor(consMEs[[set]]$data, mats, use = "p", method = cor_method)

        # Get P-value matrix
        pmat <- WGCNA::corPvalueFisher(cormat, consensus$exprSize$nSamples[set])

        results <- list(cor = cormat, pvals = pmat)
        return(results)
    })

    # Store correlation and P-value matrices in separate objects
    cormats <- lapply(matrices, function(x) return(x$cor))
    pmats <- lapply(matrices, function(x) return(x$pvals))

    # Create matrices of consensus correlation and P-values
    cons_cor <- matrix(
        NA, nrow(cormats[[1]]), ncol(cormats[[1]]),
        dimnames = list(rownames(cormats[[1]]), colnames(cormats[[1]]))
    )
    cons_pval <- cons_cor

    ## If there is a difference in sign for element m[i,j], add FALSE to it
    neg <- Reduce(`&`, lapply(cormats, function(x) return(x < 0)))
    pos <- Reduce(`&`, lapply(cormats, function(x) return(x > 0)))

    ## For every element, keep minimum cor and maximum P-value of all sets
    cons_cor[neg] <- do.call(pmax, lapply(cormats, function(x) return(x[neg])))
    cons_cor[pos] <- do.call(pmin, lapply(cormats, function(x) return(x[pos])))
    cons_pval[neg] <- do.call(pmax, lapply(pmats, function(x) return(x[neg])))
    cons_pval[pos] <- do.call(pmax, lapply(pmats, function(x) return(x[pos])))

    # Reshape to long format and merge data frames into one
    v <- c("ME", "trait")
    cor_long <- reshape2::melt(cons_cor, value.name = "cor", varnames = v)
    p_long <- reshape2::melt(cons_pval, value.name = "pvalue", varnames = v)
    final_df <- merge(cor_long, p_long)

    # Add a column `group` with metadata variable name
    var_levels <- Reduce(rbind, lapply(seq_along(metadata[[1]]), function(x) {
        return(data.frame(
            trait = unique(metadata[[1]][, x]),
            group = names(metadata[[1]])[x]
        ))
    }))
    final_df <- merge(final_df, var_levels)

    return(final_df)
}


