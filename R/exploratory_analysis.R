
#' Plot heatmap of hierarchically clustered sample correlations or gene expression
#'
#' @param exp A gene expression data frame with genes in row names
#' and samples in column names or a `SummarizedExperiment` object.
#' @param col_metadata A data frame containing sample names in row names and
#' sample annotation in the subsequent columns. The maximum number of columns
#' is 2 to ensure legends can be visualized. Ignored if `exp` is
#' a `SummarizedExperiment` object, since the function will extract colData.
#' Default: NA.
#' @param row_metadata A data frame containing gene IDs in row names and
#' gene functional classification in the first column. Only one column
#' is allowed to ensure legends can be visualized. Default: NA.
#' @param cor_method Correlation method. One of 'spearman' or 'pearson'.
#' Default is 'spearman'.
#' @param type Type of heatmap to plot. One of 'samplecor' (sample correlations)
#' or 'expr'. Default: 'samplecor'.
#' @param palette RColorBrewer palette to use. Default is "Blues" for sample
#' correlation heatmap and "YlOrRd" for gene expression heatmap.
#' @param log_trans Logical indicating whether to log transform the expression data or not. Default: FALSE.
#' @param cluster_rows Logical indicating whether to cluster rows or not.
#' Default: TRUE.
#' @param cluster_cols Logical indicating whether to cluster columns or not.
#' Default: TRUE.
#' @param show_rownames Logical indicating whether to show row names or not.
#' Default: FALSE.
#' @param show_colnames Logical indicating whether to show column names or not.
#' Default is TRUE.
#' @param scale Character indicating if values should be centered and scaled in
#' rows, columns, or none. One of 'row', 'column', or 'none'. Default: 'none'.
#' @param fontsize Base font size for the plot.
#' @param cutree_rows Number of clusters into which rows are divided.
#' Default: NA (no division).
#' @param cutree_cols Number of clusters into which columns are divided.
#' Default: NA (no division).
#' @param ... Additional arguments to be passed to \code{ComplexHeatmap::pheatmap()}.
#'
#' @return A heatmap of sample correlations or gene expression.
#' @author Fabricio Almeida-Silva
#' @seealso
#'  \code{\link[RColorBrewer]{RColorBrewer}}
#' @rdname plot_heatmap
#' @export
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ComplexHeatmap pheatmap
#' @importFrom methods is
#' @examples
#' \donttest{
#' data(filt.se)
#' plot_heatmap(filt.se)
#' }
plot_heatmap <- function(exp, col_metadata = NA, row_metadata = NA,
                         cor_method = 'spearman', type = "samplecor",
                         palette = NULL, log_trans = FALSE,
                         cluster_rows = TRUE, cluster_cols = TRUE,
                         show_rownames = FALSE, show_colnames = TRUE,
                         scale = "none", fontsize = 9,
                         cutree_rows = NA, cutree_cols = NA, ...) {
    fexp <- handleSE(exp)
    if(methods::is(exp, "SummarizedExperiment")) {
        col_metadata <- as.data.frame(SummarizedExperiment::colData(exp))
    }

    # Rename and reorder columns based on annotation and handle colors
    col_metadata <- sample_cols_heatmap(col_metadata, fexp)$col_metadata
    fexp <- sample_cols_heatmap(col_metadata, fexp)$fexp
    annotation_color <- sample_cols_heatmap(col_metadata, fexp)$annotation_colors

    # Rename and reorder rows based on annotation and handle colors
    row_metadata <- gene_cols_heatmap(row_metadata, fexp,
                                      annotation_color)$row_metadata
    fexp <- gene_cols_heatmap(row_metadata, fexp, annotation_color)$fexp
    annotation_color <- gene_cols_heatmap(row_metadata, fexp,
                                          annotation_color)$annotation_color

    if(log_trans) {
        fexp <- log2(fexp+1)
    }

    if(is.null(palette)) {
        if(type == "samplecor") {
            pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100)
        } else if(type == "expr") {
            pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(100)
        } else {
            stop("Please, specify a valid type. One of 'samplecor' or 'expr'.")
        }
    } else {
        pal <- colorRampPalette(RColorBrewer::brewer.pal(9, palette))(100)
    }


    if(type == "samplecor") {
        x <- cor(fexp, method = cor_method)
        title <- "Pairwise correlations between samples"
    } else if(type == "expr") {
        x <- fexp
        title <- "Gene expression heatmap"
    } else {
        stop("Please, specify a type. One of 'samplecor' or 'expr'.")
    }

    hm <- ComplexHeatmap::pheatmap(as.matrix(x), color=pal, border_color = NA,
                                   show_rownames = show_rownames,
                                   show_colnames = show_colnames,
                                   annotation_row = row_metadata,
                                   annotation_col = col_metadata,
                                   cluster_rows = cluster_rows,
                                   cluster_cols = cluster_cols,
                                   scale = scale, fontsize = fontsize,
                                   main=title, cutree_rows = cutree_rows,
                                   cutree_cols = cutree_cols,
                                   annotation_colors = annotation_color, ...)
    return(hm)

}

#' Plot Principal Component Analysis (PCA) of samples
#'
#' @param exp A gene expression data frame with genes in row names
#' and samples in column names or a `SummarizedExperiment` object.
#' @param metadata A data frame containing sample names in row names
#' and sample annotation in the first column. Ignored if `exp` is
#' a `SummarizedExperiment` object, since the function will extract colData.
#' @param log_trans Logical. If TRUE, the expression data frame will be log
#' transformed with log2(exp+1).
#' @param PCs Principal components to be plotted on the x-axis and y-axis,
#' respectively. One of "1x2", "1x3" or "2x3. Default is "1x2".
#' @param size Numeric indicating the point size. Default is 2.
#' @return A ggplot object with the PCA plot.
#' @author Fabricio Almeida-Silva
#' @seealso
#'  \code{\link[ggplot2]{ggplot}}
#' @rdname plot_PCA
#' @export
#' @importFrom ggplot2 ggplot aes aes_ geom_point scale_color_manual labs theme_classic ggtitle theme element_text
#' @importFrom SummarizedExperiment colData
#' @examples
#' data(zma.se)
#' plot_PCA(zma.se, log_trans = TRUE)
plot_PCA <- function(exp, metadata, log_trans = FALSE, PCs = "1x2", size = 2) {
    fexp <- handleSE(exp)
    if(is(exp, "SummarizedExperiment")) {
        metadata <- as.data.frame(SummarizedExperiment::colData(exp))
    }
    if(ncol(metadata) > 1) { stop("Sample metadata must contain one column.") }
    if(log_trans) {
        pca <- prcomp(t(log2(fexp+1)))
    } else {
        pca <- prcomp(t(fexp))
    }

    cols <- custom_palette(1)[seq_along(unique(metadata[,1]))]
    pca_df <- as.data.frame(pca$x)
    pca_df$`Sample group` <- metadata[rownames(pca_df), 1]
    pca_df$SampleID <- rownames(pca_df)
    var_explained <- as.data.frame(round(100 * pca$sdev^2 / sum(pca$sdev^2), 1))
    rownames(var_explained) <- colnames(pca$x)

    if(PCs == "1x2") {
        aes_map <- ggplot2::ggplot(pca_df, ggplot2::aes_(~PC1, ~PC2, color = ~`Sample group`))
        labs <- ggplot2::labs(x = paste("PC1 (", var_explained["PC1", ], "%)", sep = ""),
                              y = paste("PC2 (", var_explained["PC2", ], "%)", sep = ""))
    } else if(PCs == "1x3") {
        aes_map <- ggplot2::ggplot(pca_df, ggplot2::aes_(~PC1, ~PC3, color = ~`Sample group`))
        labs <- ggplot2::labs(x = paste("PC1 (", var_explained["PC1", ], "%)", sep = ""),
                              y = paste("PC3 (", var_explained["PC3", ], "%)", sep = ""))
    } else if(PCs == "2x3") {
        aes_map <- ggplot2::ggplot(pca_df, ggplot2::aes_(~PC2, ~PC3, color = ~`Sample group`))
        labs <- ggplot2::labs(x = paste("PC2 (", var_explained["PC2", ], "%)", sep = ""),
                              y = paste("PC3 (", var_explained["PC3", ], "%)", sep = ""))
    } else {
        stop("Please, specify the PCs to be plotted. One of '1x2', '1x3', or '2x3'.")
    }
    p <- aes_map +
        ggplot2::geom_point(size = size) +
        ggplot2::scale_color_manual(values = cols) +
        labs +
        ggplot2::theme_classic() +
        ggplot2::ggtitle("Principal component analysis of samples") +
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
    return(p)
}


#' Get housekeeping genes from global expression profile
#'
#' @param exp A gene expression data frame with genes in row names and
#' samples in column names or a `SummarizedExperiment` object.
#'
#' @details This function identifies housekeeping genes, which are
#' broadly expressed genes with low variation in a global scale across samples.
#' For some cases, users would want to remove these genes as they are not
#' interesting for coexpression network analyses.
#' See references for more details.
#'
#' @return Character vector of housekeeping gene IDs.
#' @author Fabricio Almeida-Silva
#' @rdname get_HK
#' @export
#' @references
#' Machado, F.B., Moharana, K.C., Almeida‐Silva, F., Gazara, R.K.,
#' Pedrosa‐Silva, F., Coelho, F.S., Grativol, C. and Venancio, T.M. (2020),
#' Systematic analysis of 1298 RNA‐Seq samples and construction of a
#' comprehensive soybean (Glycine max) expression atlas. Plant J, 103: 1894-1909.
#' @examples
#' data(zma.se)
#' hk <- get_HK(zma.se)
get_HK <- function(exp) {
    exp <- handleSE(exp)
    exp <- exp
    exp[exp < 1] <- 0 #expression values below 1 are considered as not expressed
    ncols <- ncol(exp)
    final.exp <- exp[rowSums(exp > 0) == ncols,]
    final.exp$mean <- rowMeans(final.exp[, seq_len(ncols)])
    final.exp$sd <-  apply(final.exp[, seq_len(ncols)], 1, sd)
    final.exp$covar <- final.exp$sd/final.exp$mean
    final.exp$max <- apply(final.exp[, seq_len(ncols)], 1, max)
    final.exp$min <- apply(final.exp[, seq_len(ncols)], 1, min)
    final.exp$MFC <- final.exp$max/final.exp$min
    final.exp$MFC.CoV <- final.exp$MFC * final.exp$covar
    hk <- head(final.exp[order(final.exp$MFC.CoV, decreasing=FALSE),],
               n = nrow(final.exp)*0.25)
    hk.genes <- rownames(hk)
    return(hk.genes)
}

#' Plot expression profile of given genes across samples
#'
#' @param genes Character vector containing a subset of genes from which
#' edges will be extracted. It can be ignored if \code{plot_module} is TRUE.
#' @param exp A gene expression data frame with genes in row names and
#' samples in column names or a `SummarizedExperiment` object.
#' @param metadata A data frame containing sample names in row names
#' and sample annotation in the first column. Ignored if `exp`
#' is a `SummarizedExperiment` object, since the function will extract colData.
#' @param plot_module Logical indicating whether to plot a whole module or not.
#' If set to FALSE, \code{genes} must be specified.
#' @param net List object returned by \code{exp2net}.
#' @param modulename Name of the module to plot.
#'
#' @return A ggplot object showing the expression profile of some genes
#' across all samples.
#' @author Fabricio Almeida-Silva
#' @seealso
#'  \code{\link[reshape2]{melt}}
#' @rdname plot_expression_profile
#' @export
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes_ geom_tile aes geom_line stat_summary ggtitle
#' @importFrom SummarizedExperiment colData
#' @examples
#' data(zma.se)
#' data(filt.se)
#' genes <- rownames(filt.se)
#' plot_expression_profile(genes=genes, exp=zma.se, plot_module=FALSE)
plot_expression_profile <- function(genes, exp, metadata, plot_module = TRUE,
                                    net, modulename) {
    if(is(exp, "SummarizedExperiment")) {
        metadata <- as.data.frame(SummarizedExperiment::colData(exp))
    }
    exp <- handleSE(exp)

    title <- "Expression profile"
    if(plot_module) {
        genes_modules <- net$genes_and_modules
        genes <- genes_modules[genes_modules$Modules == modulename, "Genes"]
        title <- paste("Expression profile for module", modulename)
    }
    fexp <- as.data.frame(exp[rownames(exp) %in% genes, ])
    fexp$id <- rownames(fexp)
    melt_exp <- reshape2::melt(fexp, "id", variable.name = "Samples",
                               value.name = "Expression")

    # Reorder samples by tissue
    metadata <- metadata[colnames(fexp)[-ncol(fexp)], , drop=FALSE]
    metadata$Sample <- rownames(metadata)
    colnames(metadata) <- c("Sample group", "Sample")
    metadata <- metadata[order(metadata$`Sample group`), ]
    metadata$Sample <- factor(metadata$Sample, levels = metadata$Sample)
    melt_exp[, "Samples"] <- factor(melt_exp[, "Samples"], levels = metadata$Sample)

    # Background tiles
    background <- mean(melt_exp[, "Expression"])

    # Plot expression profiles
    cols <- custom_palette(1)[seq_along(unique(metadata$`Sample group`))]
    p <- ggplot2::ggplot(melt_exp, ggplot2::aes_(x = ~Samples, y = ~Expression)) +
        ggplot2::geom_tile(data = metadata, alpha = 0.4, height = Inf,
                           ggplot2::aes_(x = ~Sample, y = ~background,
                      fill = ~`Sample group`)) +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::geom_line(ggplot2::aes_(group = ~id), alpha = 0.2,
                       color = "firebrick") +
        ggplot2::stat_summary(ggplot2::aes(group = 1), size = 1,
                              fun = "median", geom = "line") +
        theme_exp_profile() +
        ggplot2::ggtitle(title)

    return(p)
}


#' Plot number of genes per module
#'
#' @param net List object returned by \code{exp2gcn}.
#' @return A ggplot object with a bar plot of gene number in each module.
#' @seealso
#'  \code{\link[ggpubr]{ggbarplot}}
#' @rdname plot_ngenes_per_module
#' @export
#' @importFrom ggpubr ggbarplot
#' @importFrom ggplot2 theme element_text
#' @examples
#' data(filt.se)
#' gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson")
#' plot_ngenes_per_module(gcn)
plot_ngenes_per_module <- function(net = NULL) {
    genes_and_modules <- net$genes_and_modules
    frequency_df <- as.data.frame(table(genes_and_modules$Modules),
                                  stringsAsFactors=FALSE)
    names(frequency_df) <- c("Module", "Frequency")
    frequency_df <- frequency_df[order(frequency_df$Frequency, decreasing = TRUE), ]
    cols <- unique(frequency_df$Module)
    ggpubr::ggbarplot(data=frequency_df, x="Module", y="Frequency",
                      fill="Module", palette=cols,
                      legend="none", title="Number of genes per module",
                      x.text.angle=60, font.title="bold",
                      label=TRUE, label.pos="out", ylim=c(0, NA)) +
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
}


