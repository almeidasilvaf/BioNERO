

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
#' @importFrom ggplot2 ggplot aes_ geom_point labs theme_bw scale_color_manual
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
    p1 <- ggplot(sft_df, aes_(x = ~power, y = ~fit)) +
        geom_point(aes_(color = ~Set)) +
        ggrepel::geom_text_repel(aes_(label = ~power, color = ~Set)) +
        scale_color_manual(values = cols) +
        labs(
            x = "Soft threshold (power)",
            y = expression(paste("Scale-free topology fit - ", R^{2})),
            title = "Scale independence"
        ) +
        theme_bw() +
        ylim(c(0,1)) +
        geom_hline(yintercept = rsquared, color="brown3") +
        theme(legend.position = "none")

    # Plot 2
    p2 <- ggplot(sft_df, aes_(x = ~power, y = ~meank)) +
        geom_point(aes_(color = ~Set)) +
        ggrepel::geom_text_repel(aes_(label = ~power, color = ~Set),
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
#' @param verbose Logical indicating whether to display progress
#' messages or not. Default: FALSE.
#'
#' @return A list containing 4 elements: \describe{
#'   \item{consModules}{Consensus module assignments}
#'   \item{consMEs}{Consensus module eigengenes}
#'   \item{exprSize}{Description of the multi-set object returned by the function \code{WGCNA::checkSets}}
#'   \item{sampleInfo}{Metadata for each expression set}
#'   \item{genes_cmodules}{Data frame of genes and consensus modules}
#'   \item{dendro_plot_objects}{Objects to be used in dendrogram plotting}
#' }
#'
#' @seealso
#'  \code{\link[dynamicTreeCut]{cutreeDynamic}}
#' @rdname consensus_modules
#' @export
#' @importFrom WGCNA checkSets adjacency TOMsimilarity pquantile labels2colors multiSetMEs consensusMEDissimilarity mergeCloseModules
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
consensus_modules <- function(exp_list, metadata, power, cor_method = "spearman",
                              net_type = "signed hybrid",
                              module_merging_threshold = 0.8,
                              verbose = FALSE) {
    metadata <- handle_metadata(exp_list, metadata)
    exp_list <- handleSElist(exp_list)
    nSets <- length(exp_list)

    # Build multi-set object
    multiExp <- lapply(exp_list, function(x) {
        element <- list(data=as.data.frame(t(x)))
        return(element)
    })

    # Check if the data has the correct format for downstream analysis
    expSize <- WGCNA::checkSets(multiExp)
    nGenes <- expSize$nGenes
    nSamples <- expSize$nSamples

    # Build multi-set object with metadata
    sampleinfo <- lapply(seq_along(multiExp), function(x) {
        sinfo <- metadata[rownames(multiExp[[x]]$data), , drop=FALSE]
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
    TOM <- lapply(adj, function(x) {
        if(net_type == "signed hybrid") {
            tom <- WGCNA::TOMsimilarity(x, TOMType = "signed")
        } else if(net_type == "signed") {
            tom <- WGCNA::TOMsimilarity(x, TOMType = "signed Nowick")
        } else {
            tom <- WGCNA::TOMsimilarity(x, TOMType = "unsigned")
        }
        return(tom)
    })

    # Scaling TOM to make them comparable
    scaleP <- 0.95
    nSamples <- as.integer(1 / (1-scaleP) * 100)
    scaleSample <- sample(nGenes * (nGenes-1) / 2, size = nSamples)

    scaledTOM <- lapply(seq_len(nSets), function(x) {
        TOMScalingSamples <- as.dist(TOM[[x]])[scaleSample]
        scaleQuant <- quantile(TOMScalingSamples, probs = scaleP, type = 8)
        tom_scaled <- TOM[[x]]
        if(x > 1) {
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
        distM = 1-consensusTOM,
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
    consMEDiss <- WGCNA::consensusMEDissimilarity(unmergedMEs)
    consMETree <- hclust(as.dist(consMEDiss), method = "average")

    merging_threshold <- 1-module_merging_threshold
    merge <- WGCNA::mergeCloseModules(
        multiExp, unmergedColors, cutHeight = merging_threshold, verbose = 0
    )

    moduleLabels <- merge$colors
    moduleColors <- WGCNA::labels2colors(moduleLabels, colorSeq = palette)
    genes_cmod <- data.frame(
        Genes = rownames(adj[[1]]), Cons_modules = moduleColors
    )
    consMEs <- merge$newMEs

    # Plot dendrogram with merged colors
    result_list <- list(
        consModules = moduleColors, consMEs = consMEs, exprSize = expSize,
        sampleInfo = sampleinfo, genes_cmodules = genes_cmod,
        dendro_plot_objects = list(tree = consTree, unmerged = unmergedColors)
    )
    return(result_list)
}

#' Plot dendrogram of genes and consensus modules
#'
#' @param consensus Consensus network returned by \code{consensus_modules}.
#'
#' @return A base plot with the gene dendrogram and modules.
#' @importFrom WGCNA plotDendroAndColors
#' @importFrom graphics par layout
#' @export
#' @rdname plot_dendro_and_cons_colors
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
#' plot_dendro_and_cons_colors(cons_mod)
plot_dendro_and_cons_colors <- function(consensus) {
    on.exit(graphics::layout(1))
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar), add = TRUE, after = FALSE)
    WGCNA::plotDendroAndColors(
        consensus$dendro_plot_objects$tree,
        cbind(
            consensus$dendro_plot_objects$unmerged,
            consensus$consModules
        ),
        c("Unmerged", "Merged"),
        dendroLabels = FALSE, hang = 0.03,
        addGuide = TRUE, guideHang = 0.05)
}


#' Correlate set-specific modules and consensus modules to sample information
#'
#' @param consensus Consensus network returned by \code{consensus_modules}.
#' @param cor_method Correlation method to be used. One of 'spearman' or
#' 'pearson'. Default is 'spearman'.
#' @param continuous_trait Logical indicating if trait is a continuous variable.
#' Default is FALSE.
#' @param palette RColorBrewer's color palette to use. Default is "RdYlBu",
#' a palette ranging from blue to red.
#' @param cex.lab.x Font size for x axis labels. Default: 0.6.
#' @param cex.lab.y Font size for y axis labels. Default: 0.6.
#' @param cex.text Font size for numbers inside matrix. Default: 0.6.
#' @param transpose Logical indicating whether to transpose the heatmap of not.
#' Default is FALSE.
#'
#' @return Data frame of consensus module-trait correlations and p-values.
#' @details Significance levels:
#' 1 asterisk: significant at alpha = 0.05.
#' 2 asterisks: significant at alpha = 0.01.
#' 3 asterisks: significant at alpha = 0.001.
#' no asterisk: not significant.
#'
#' @seealso
#'  \code{\link[WGCNA]{corPvalueFisher}},\code{\link[WGCNA]{labels2colors}},\code{\link[WGCNA]{labeledHeatmap}},\code{\link[WGCNA]{blueWhiteRed}}
#' @rdname consensus_trait_cor
#' @export
#' @importFrom WGCNA corPvalueFisher labels2colors labeledHeatmap
#' @importFrom reshape2 melt
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics par layout
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
consensus_trait_cor <- function(consensus, cor_method = "spearman",
                                continuous_trait = FALSE,
                                palette="RdYlBu",
                                cex.lab.x=0.6, cex.lab.y=0.6,
                                cex.text=0.6,
                                transpose=FALSE) {

    consMEs <- consensus$consMEs
    exprSize <- consensus$exprSize
    sampleInfo <- consensus$sampleInfo
    sampleInfo <- lapply(
        sampleInfo, handle_trait_type, continuous_trait = continuous_trait
    )
    nSets <- exprSize$nSets

    # Calculate the correlations
    mod_trait_cor <- lapply(seq_len(nSets), function(set) {
        if(cor_method == "spearman") {
            moduleTraitCor <- cor(
                consMEs[[set]]$data, sampleInfo[[set]], use = "p",
                method = "spearman"
            )
            moduleTraitPvalue <- WGCNA::corPvalueFisher(
                moduleTraitCor, exprSize$nSamples[set]
            )
        } else if(cor_method == "pearson") {
            moduleTraitCor <- cor(
                consMEs[[set]]$data, sampleInfo[[set]], use = "p",
                method = "pearson"
            )
            moduleTraitPvalue <- WGCNA::corPvalueFisher(
                moduleTraitCor, exprSize$nSamples[set]
            )
        }
        results <- list(moduleTraitCor, moduleTraitPvalue)
        return(results)
    })
    moduleTraitCor <- lapply(mod_trait_cor, function(x) return(x[[1]]))
    moduleTraitPvalue <- lapply(mod_trait_cor, function(x) return(x[[2]]))

    MEColors <- substring(names(consMEs[[1]]$data), 3)
    MEColorNames <- paste("ME", MEColors, sep="")

    # Plot the module-trait relationship table for sets
    cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(10, palette)))(100)

    # Initialize matrices to hold the consensus correlation and p-value
    cons_cor <- matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
    cons_pval <- matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))

    neg <- lapply(moduleTraitCor, function(x) return(x < 0))
    neg <- Reduce(`&`, neg)
    pos <- lapply(moduleTraitCor, function(x) return(x > 0))
    pos <- Reduce(`&`, pos)

    cons_cor[neg] <- do.call(pmax, lapply(moduleTraitCor, function(x) return(x[neg])))
    cons_cor[pos] <- do.call(pmin, lapply(moduleTraitCor, function(x) return(x[pos])))
    cons_pval[neg] <- do.call(pmax, lapply(moduleTraitPvalue, function(x) return(x[neg])))
    cons_pval[pos] <- do.call(pmax, lapply(moduleTraitPvalue, function(x) return(x[pos])))

    # Create data frame of correlations and p-values to return
    colnames(cons_cor) <- colnames(moduleTraitCor[[1]])
    rownames(cons_cor) <- rownames(moduleTraitCor[[1]])
    colnames(cons_pval) <- colnames(moduleTraitCor[[1]])
    rownames(cons_pval) <- rownames(moduleTraitCor[[1]])
    cor_long <- reshape2::melt(cons_cor)
    pval_long <- reshape2::melt(cons_pval)
    combined_long <- merge(cor_long, pval_long, by = c("Var1", "Var2"))
    colnames(combined_long) <- c("ME", "trait", "cor", "pvalue")
    combined_long$ME <- as.character(combined_long$ME)
    combined_long$trait <- as.character(combined_long$trait)

    set <- 1
    setLabels <- c("Set1", "Set2")
    modtraitsymbol <- pval2symbol(cons_pval)
    textMatrix <- paste(signif(cons_cor, 2), modtraitsymbol, sep = "")
    textMatrix[textMatrix == "NANA"] <- "-"
    dim(textMatrix) <- dim(moduleTraitCor[[set]])
    if(transpose) {
        cons_cor <- t(cons_cor)
        textMatrix <- t(textMatrix)
        yLabels <- names(sampleInfo[[set]])
        xLabels <- MEColorNames
        xSymbols <- MEColorNames
        ySymbols <- NULL
        xColorLabels <- TRUE
        par(mar = c(5, 5, 1, 1))
    } else {
        par(mar = c(6, 8.5, 3, 3))
        xLabels <- names(sampleInfo[[set]])
        yLabels <- MEColorNames
        ySymbols <- MEColorNames
        xColorLabels <- FALSE
        xSymbols <- NULL
        xColorLabels <- FALSE
    }
    on.exit(graphics::layout(1))
    opar <- par(no.readonly=TRUE)
    on.exit(par(opar), add=TRUE, after=FALSE)
    hm <- WGCNA::labeledHeatmap(
        Matrix = cons_cor, yLabels = yLabels, xLabels = xLabels,
        ySymbols = ySymbols, xSymbols = xSymbols,
        colorLabels = FALSE, colors = cols,
        textMatrix = textMatrix, setStdMargins = FALSE,
        cex.text = cex.text, cex.lab.x = cex.lab.x, cex.lab.y = cex.lab.y,
        zlim = c(-1,1), cex.main = 1,
        main = "Consensus module-trait relationships"
    )
    return(combined_long)
}


