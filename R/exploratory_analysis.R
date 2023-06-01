
#' Plot heatmap of hierarchically clustered sample correlations or gene expression
#'
#' @param exp A gene expression data frame with genes in row names
#' and samples in column names or a `SummarizedExperiment` object.
#' @param col_metadata A data frame containing sample names in row names and
#' sample annotation in the subsequent columns. The maximum number of columns
#' is 3 to ensure legends can be visualized. Ignored if `exp` is
#' a `SummarizedExperiment` object, since the function will extract colData.
#' Default: NA.
#' @param row_metadata A data frame containing gene IDs in row names and
#' gene functional classification in the first column. The maximum number
#' of columns is 3 to ensure legends can be visualized. Default: NA.
#' @param coldata_cols A vector (either numeric or character) indicating
#' which columns should be extracted from column metadata if \strong{exp}
#' is a `SummarizedExperiment` object. The vector can contain column
#' indices (numeric) or column names (character). By default, all columns are
#' used.
#' @param rowdata_cols A vector (either numeric or character) indicating
#' which columns should be extracted from row metadata if \strong{exp}
#' is a `SummarizedExperiment` object. The vector can contain column
#' indices (numeric) or column names (character). By default, all columns are
#' used.
#' @param type Type of heatmap to plot. One of 'samplecor' (sample correlations)
#' or 'expr'. Default: 'samplecor'.
#' @param cor_method Correlation method to use in
#' case \strong{type} is "samplecor". One of 'spearman' or 'pearson'.
#' Default is 'spearman'.
#' @param palette RColorBrewer palette to use. Default is "Blues" for sample
#' correlation heatmaps and "YlOrRd" for gene expression heatmaps.
#' @param log_trans Logical indicating whether to log transform the expression
#' data or not. Default: FALSE.
#' @param ... Additional arguments to be passed
#' to \code{ComplexHeatmap::pheatmap()}. These arguments can be used to control
#' heatmap aesthetics, such as show/hide row and column names,
#' change font size, activate/deactivate hierarchical clustering, etc. For a
#' complete list of the options, see \code{?ComplexHeatmap::pheatmap()}.
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
                         coldata_cols = NULL, rowdata_cols = NULL,
                         type = "samplecor", cor_method = 'spearman',
                         palette = NULL, log_trans = FALSE, ...) {

    fexp <- handleSE(exp)
    if(methods::is(exp, "SummarizedExperiment")) {
        col_metadata <- se2metadata(exp, coldata_cols = coldata_cols)$coldata
        row_metadata <- se2metadata(exp, rowdata_cols = rowdata_cols)$rowdata
    }

    # Get column metadata with sorted rows + list of named vectors with colors
    cdata <- metadata2colors(col_metadata)
    col_metadata <- cdata$metadata
    col_colors <- cdata$colors

    # Get row metadata with sorted rows + list of named vectors with colors
    rdata <- metadata2colors(row_metadata)
    row_metadata <- rdata$metadata
    row_colors <- rdata$colors

    annotation_colors <- c(col_colors, row_colors)

    # Reorder rows and columns of the expression matrix based on metadata
    if(is.data.frame(col_metadata)) { fexp <- fexp[, rownames(col_metadata)] }
    if(is.data.frame(row_metadata)) { fexp <- fexp[rownames(row_metadata), ] }

    # Should data be log-transformed?
    if(log_trans) { fexp <- log2(fexp + 1) }

    # Get attributes (color palette, matrix, title, and legend title)
    hm_attr <- heatmap_attributes(fexp, palette, type, cor_method)

    # Plot heatmap
    hm <- ComplexHeatmap::pheatmap(
        as.matrix(hm_attr$mat),
        name = hm_attr$name,
        color = hm_attr$pal,
        border_color = NA,
        annotation_row = row_metadata,
        annotation_col = col_metadata,
        main = hm_attr$title,
        annotation_colors = annotation_colors,
        ...
    )

    return(hm)
}

#' Plot Principal Component Analysis (PCA) of samples
#'
#' @param exp A gene expression data frame with genes in row names
#' and samples in column names or a `SummarizedExperiment` object.
#' @param metadata A data frame of sample metadata containing sample names
#' in row names and sample annotation in subsequent columns.
#' Ignored if `exp` is a `SummarizedExperiment` object, since colData will be
#' automatically extracted.
#' @param metadata_cols A vector (either numeric or character) indicating
#' which columns should be extracted from column metadata if \strong{exp}
#' is a `SummarizedExperiment` object. The vector can contain column
#' indices (numeric) or column names (character). By default, all columns are
#' used.
#' @param log_trans Logical indicating whether the gene expression matrix
#' should be log transformed using \code{log(exp + 1)}. Default: FALSE.
#' @param PCs Numeric vector of length 2 indicating the principal components
#' to be plotted on the x-axis and y-axis, respectively.
#' Default: \code{c(1, 2)}.
#' @param size Numeric indicating the point size. Default is 2.
#' @return A ggplot object with the PCA plot.
#'
#' @author Fabricio Almeida-Silva
#' @seealso
#'  \code{\link[ggplot2]{ggplot}}
#' @rdname plot_PCA
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs
#' theme_classic ggtitle theme element_text
#' @importFrom SummarizedExperiment colData
#' @importFrom stats prcomp
#' @importFrom rlang .data
#' @examples
#' data(zma.se)
#' plot_PCA(zma.se, log_trans = TRUE)
plot_PCA <- function(exp, metadata, metadata_cols = NULL,
                     log_trans = FALSE, PCs = c(1, 2), size = 2) {

    # Get expression matrix and colData if exp is a `SummarizedExperiment`
    fexp <- handleSE(exp)
    if(methods::is(exp, "SummarizedExperiment")) {
        metadata <- se2metadata(exp, coldata_cols = metadata_cols)$coldata
    }

    if(ncol(metadata) > 2) { stop("Sample metadata must contain up to 2 columns.") }
    if(log_trans) { fexp <- log2(fexp + 1) }

    # Perform PCA
    pca <- prcomp(t(fexp))

    # Get 2 data frames: PCs with metadata, and variance explained
    pcs <- merge(as.data.frame(pca$x), metadata, by = "row.names")
    varexp <- data.frame(
        row.names = colnames(pca$x),
        Var = round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
    )

    # Get plot parameters
    pc_x <- paste0("PC", PCs[1])
    pc_y <- paste0("PC", PCs[2])
    colors <- custom_palette(1)[seq_along(unique(metadata[, 1]))]
    color_by <- names(metadata)[1]
    shape_by <- names(metadata)[2]

    point <- geom_point(aes(color = .data[[color_by]]), size = size)
    if(ncol(metadata) > 1) {
        point <- geom_point(
            aes(color = .data[[color_by]], shape = .data[[shape_by]]),
            size = size
        )
    }

    # Create plot
    p <- ggplot(pcs, aes(x = .data[[pc_x]], y = .data[[pc_y]])) +
        point +
        labs(
            title = "PCA of samples",
            x = paste0(pc_x, " (", varexp[pc_x, ], "%)"),
            y = paste0(pc_y, " (", varexp[pc_y, ], "%)")
        ) +
        scale_color_manual(values = colors) +
        theme_classic()

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
    exp[exp < 1] <- 0 # expression values < 1 are considered as not expressed
    ncols <- ncol(exp)
    final.exp <- exp[rowSums(exp > 0) == ncols,]
    final.exp$mean <- rowMeans(final.exp[, seq_len(ncols)])
    final.exp$sd <-  apply(final.exp[, seq_len(ncols)], 1, sd)
    final.exp$covar <- final.exp$sd / final.exp$mean
    final.exp$max <- apply(final.exp[, seq_len(ncols)], 1, max)
    final.exp$min <- apply(final.exp[, seq_len(ncols)], 1, min)
    final.exp$MFC <- final.exp$max / final.exp$min
    final.exp$MFC.CoV <- final.exp$MFC * final.exp$covar
    hk <- head(
        final.exp[order(final.exp$MFC.CoV, decreasing = FALSE),],
        n = nrow(final.exp) * 0.25
    )
    hk.genes <- rownames(hk)
    return(hk.genes)
}

#' Plot expression profile of given genes across samples
#'
#' @param genes Character vector containing a subset of genes from which
#' edges will be extracted. It can be ignored if \code{plot_module} is TRUE.
#' @param exp A gene expression data frame with genes in row names and
#' samples in column names or a `SummarizedExperiment` object.
#' @param metadata A data frame of sample metadata containing sample names
#' in row names and sample annotation in subsequent columns.
#' Ignored if `exp` is a `SummarizedExperiment` object, since colData will be
#' automatically extracted.
#' @param metadata_cols A character or numeric scalar indicating
#' which column should be extracted from column metadata if \strong{exp}
#' is a `SummarizedExperiment` object. The column to be extracted can be
#' represented by indices (numeric) or column names (character).
#' By default, the first column is used.
#' @param plot_module Logical indicating whether to plot a whole module or not.
#' If set to FALSE, \code{genes} must be specified.
#' @param net List object returned by \code{exp2gcn}.
#' @param modulename Name of the module to plot.
#' @param bg_line Character indicating what to show in the background (black)
#' line. One of "mean" or "median". Default: "mean".
#'
#' @return A ggplot object showing the expression profile of some genes
#' across all samples.
#' @author Fabricio Almeida-Silva
#'
#' @rdname plot_expression_profile
#' @export
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_tile aes geom_line stat_summary labs
#' theme theme_classic element_blank
#' @examples
#' data(zma.se)
#' data(filt.se)
#' genes <- rownames(filt.se)
#' plot_expression_profile(genes = genes, exp = zma.se, plot_module = FALSE)
plot_expression_profile <- function(genes, exp, metadata, metadata_cols = 1,
                                    plot_module = TRUE, net, modulename,
                                    bg_line = "mean") {

    # Extract expression matrix and colData if exp is a `SummarizedExperiment`
    if(is(exp, "SummarizedExperiment")) {
        metadata <- se2metadata(exp, coldata_cols = metadata_cols)$coldata
    }
    exp <- as.matrix(handleSE(exp))

    # Reorder rows of the metadata based on levels of the chosen variable
    cols <- metadata2colors(metadata)
    metadata <- cols$metadata
    colors <- cols$colors[[1]]

    # Get genes to plot and plotting parameters
    title <- "Expression profile"
    if(plot_module) {
        genes_modules <- net$genes_and_modules
        genes <- genes_modules[genes_modules$Modules == modulename, "Genes"]
        title <- paste("Expression profile for module", modulename)
    }

    # Reshape expression matrix to long format
    exp_long <- reshape2::melt(exp[rownames(exp) %in% genes, ])
    names(exp_long) <- c("Gene", "Sample", "Expression")

    # Convert `Sample` column to factor and reorder levels based on metadata
    metadata$Sample <- factor(rownames(metadata), levels = rownames(metadata))
    exp_long$Sample <- factor(exp_long$Sample, levels = rownames(metadata))

    # Y-axis limit for the tile
    metadata$background <- mean(exp_long$Expression)

    # Plot expression profiles
    var <- names(metadata)[1]
    p <- ggplot(exp_long, aes(x = .data$Sample, y = .data$Expression)) +
        geom_tile(data = metadata, aes(
            x = .data$Sample, y = .data$background, fill = .data[[var]]
        ), alpha = 0.3, height = Inf) +
        scale_fill_manual(values = colors) +
        geom_line(aes(group = .data$Gene), alpha = 0.3, color = "firebrick") +
        stat_summary(aes(group = 1), linewidth = 1, fun = bg_line, geom = "line") +
        labs(x = "", title = title) +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            legend.position = "bottom"
        )

    return(p)
}


#' Plot number of genes per module
#'
#' @param net List object returned by \code{exp2gcn}.
#' @return A ggplot object with a bar plot of gene number in each module.
#'
#' @rdname plot_ngenes_per_module
#' @export
#' @importFrom ggplot2 ggplot geom_col scale_fill_manual theme_bw labs
#' theme element_text geom_text xlim
#' @examples
#' data(filt.se)
#' gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson")
#' plot_ngenes_per_module(gcn)
plot_ngenes_per_module <- function(net = NULL) {

    # Create a data frame of gene counts per module
    genes_and_modules <- net$genes_and_modules
    frequency_df <- as.data.frame(table(genes_and_modules$Modules))
    names(frequency_df) <- c("Module", "Frequency")

    # Reorder `Frequency` column in descending order
    frequency_df <- frequency_df[order(frequency_df$Frequency, decreasing = TRUE), ]
    frequency_df$Module <- factor(
        frequency_df$Module,
        levels = frequency_df$Module
    )

    # Get colors and axis limits
    cols <- levels(frequency_df$Module)
    xmax <- ceiling(max(frequency_df$Frequency) / 100) * 100

    # Create plot
    p <- ggplot(frequency_df, aes(x = .data$Frequency, y = .data$Module)) +
        geom_col(aes(fill = .data$Module), color = "black") +
        geom_text(aes(label = .data$Frequency), hjust = -0.3) +
        scale_fill_manual(values = cols) +
        xlim(0, xmax) +
        theme_bw() +
        labs(title = "Number of genes per module") +
        theme(legend.position = "none")

    return(p)
}


