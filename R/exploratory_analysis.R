
#' Plot heatmap of hierarchically clustered sample correlations or gene expression
#'
#' @param exp Expression data frame with genes in rows and samples in columns.
#' @param col_metadata A data frame containing sample names in the first column and sample description in the subsequent columns. If there are more than 2 columns, additional columns will be interpreted as multiple column annotations.
#' @param row_metadata A data frame containing gene IDs in the first column and gene functional classification in the subsequent columns. If there are more than 2 columns, additional columns will be interpreted as multiple gene annotation classifications.
#' @param cor_method Correlation method. One of 'spearman' or 'pearson'. Default is 'spearman'.
#' @param type Type of heatmap to plot. One of 'samplecor' (sample correlations) or 'expr'. Default is 'samplecor'.
#' @param palette RColorBrewer palette to use. Default is "Blues" for sample correlation heatmap and "YlOrRd" for gene expression heatmap.
#' @param log_trans Logical. It specifies whether to log transform data or not. Default is FALSE.
#' @param cluster_rows Logical indicating whether to cluster rows or not. Default is TRUE.
#' @param cluster_cols Logical indicating whether to cluster columns or not. Default is TRUE.
#' @param show_rownames Logical indicating whether to show row names or not. Default is FALSE.
#' @param show_colnames Logical indicating whether to show column names or not. Default is TRUE.
#' @param scale Character indicating if values should be centered and scaled in rows, columns, or none. One of 'row', 'column', or 'none'. Default is 'none'.
#' @param fontsize Base fontsize for the plot.
#' @param cutree_rows Number of clusters into which rows are divided. Default is 1, which is equivalent to no division.
#' @param cutree_cols Number of clusters into which columns are divided. Default is 1, which is equivalent to no division.
#'
#' @return heatmap of hierarchically clustered samples with metadata information (optional)
#' @author Fabricio Almeida-Silva
#' @seealso
#'  \code{\link[RColorBrewer]{RColorBrewer}}
#'  \code{\link[pheatmap]{pheatmap}}
#' @rdname plot_heatmap
#' @export
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pheatmap pheatmap
plot_heatmap <- function(exp, col_metadata = NULL, row_metadata = NULL,
                         cor_method = 'spearman', type = "samplecor",
                         palette = NULL, log_trans = FALSE,
                         cluster_rows = TRUE, cluster_cols = TRUE,
                         show_rownames = FALSE, show_colnames = TRUE,
                         scale = "none", fontsize = 10,
                         cutree_rows = 1, cutree_cols = 1, ...) {

    if(is.null(col_metadata)) {
        col_metadata <- NULL
    } else if(ncol(col_metadata) == 2) {
        names(col_metadata) <- c("SampleID", "Sample class")
    } else if(ncol(col_metadata) == 3) {
        names(col_metadata) <- c("SampleID", "Sample class 1", "Sample class 2")
    } else if(ncol(col_metadata == 4)) {
        names(col_metadata) <- c("SampleID", "Sample class 1", "Sample class 2",
                                 "Sample class 3")
    } else {
        col_metadata <- NULL
    }

    if(log_trans == TRUE) {
        exp <- log2(exp+1)
    } else {
        exp <- exp
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
        pal <- RColorBrewer::brewer.pal(9, palette)
    }

    # How should we handle column annotation?
    if(is.null(col_metadata)) {
        annotation_col <- NULL
    } else {
        annotation_col <- col_metadata
        rownames(annotation_col) <- annotation_col[,1]
        annotation_col[,1] <- NULL
    }

    # How should we handle row annotation?
    if(is.null(row_metadata)) {
        annotation_row <- NULL
    } else {
        annotation_row <- row_metadata
        rownames(annotation_row) <- annotation_row[,1]
        annotation_row[,1] <- NULL
    }


    if(type == "samplecor") {
        x <- cor(exp, method = cor_method)
    } else if(type == "expr") {
        x <- exp
    } else {
        stop("Please, specify a type. One of 'samplecor' or 'expr'.")
    }

    title <- ifelse(type == "samplecor",
                    "Pairwise correlations between samples",
                    "Gene expression heatmap")

    # Define colors
    ccols <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF",
               "#9467BDFF", "#8C564BFF", "#E377C2FF", "#7F7F7FFF",
               "#BCBD22FF", "#17BECFFF", "#AEC7E8FF", "#FFBB78FF",
               "#98DF8AFF", "#FF9896FF", "#C5B0D5FF", "#C49C94FF",
               "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF")
    nlevels <- length(unique(as.character(annotation_col[,1])))
    if(nlevels <= 20) {
        custom_cols <- ccols[1:nlevels]
    } else {
        custom_cols <- colorRampPalette(ccols)(nlevels)
    }
    annotation_color <- list(`Sample class` = custom_cols)
    names(annotation_color$`Sample class`) <- unique(annotation_col[,1])


    pheatmap::pheatmap(x, color=pal, border_color = NA, height = 20,
                       show_rownames = show_rownames, show_colnames = show_colnames,
                       annotation_row = annotation_row, annotation_col = annotation_col,
                       cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                       scale = scale, fontsize = fontsize, main=title,
                       cutree_rows = cutree_rows, cutree_cols = cutree_cols,
                       annotation_colors = annotation_color, ...)

}

#' Plot Principal Component Analysis (PCA) of samples
#'
#' @param exp Expression data frame with genes in rownames and samples in column names.
#' @param metadata Data frame containing sample IDs in first column and sample description (e.g. tissue or treatment) in the second column.
#' @param log_trans Logical. If TRUE, the expression data frame will be log transformed by log2(exp+1).
#' @param PCs Principal Components to be plotted on the x-axis and y-axis, respectively. One of "1x2", "1x3" or "2x3. Default is "1x2".
#' @param interactive Logical indicating whether PCA should be plotted using Shiny's interactivity on RStudio. Default is FALSE.
#' @param size Numeric indicating the point size. Default is 2.
#'
#' @return A ggplot object with the PCA plot.
#' @author Fabricio Almeida-Silva
#' @seealso
#'  \code{\link[ggplot2]{ggplot}}
#' @rdname plot_PCA
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs theme_classic ggtitle theme element_text
#' @importFrom ggvis ggvis layer_points input_slider add_tooltip
plot_PCA <- function(exp, metadata, log_trans = FALSE, PCs = "1x2", size = 2, interactive = FALSE) {
    if (log_trans == FALSE) {
        pca <- prcomp(t(exp))
    } else {
        pca <- prcomp(t(log2(exp+1)))
    }

    # Define colors
    ccols <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF",
                     "#9467BDFF", "#8C564BFF", "#E377C2FF", "#7F7F7FFF",
                     "#BCBD22FF", "#17BECFFF", "#AEC7E8FF", "#FFBB78FF",
                     "#98DF8AFF", "#FF9896FF", "#C5B0D5FF", "#C49C94FF",
                     "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF")
    nlevels <- length(unique(as.character(metadata[,2])))
    if(nlevels <= 20) {
        custom_cols <- ccols
    } else {
        custom_cols <- colorRampPalette(ccols)(nlevels)
    }

    pca_df <- as.data.frame(pca$x)
    pca_df$Trait <- metadata[metadata[,1] %in% rownames(pca_df), 2]
    pca_df$SampleID <- rownames(pca_df)
    var_explained <- as.data.frame(round(100 * pca$sdev^2 / sum(pca$sdev^2), 1))
    rownames(var_explained) <- colnames(pca$x)
    if(interactive == FALSE) {
        if(PCs == "1x2") {
            p <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, color = Trait)) +
                ggplot2::geom_point(size = size) +
                ggplot2::scale_color_manual(values=custom_cols) +
                ggplot2::labs(x = paste("PC1 (", var_explained["PC1", ], "%)", sep = ""), y = paste("PC2 (", var_explained["PC2", ], "%)", sep = "")) +
                ggplot2::theme_classic() +
                ggplot2::ggtitle("Principal component analysis of samples") +
                ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
        } else if (PCs == "1x3") {
            p <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC3, color = Trait)) +
                ggplot2::geom_point(size = size) +
                ggplot2::scale_color_manual(values=custom_cols) +
                ggplot2::labs(x = paste("PC1 (", var_explained["PC1", ], "%)", sep = ""), y = paste("PC3 (", var_explained["PC3", ], "%)", sep = "")) +
                ggplot2::theme_classic() +
                ggplot2::ggtitle("Principal component analysis of samples") +
                ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
        } else if (PCs == "2x3") {
            p <- ggplot2::ggplot(pca_df, ggplot2::aes(PC2, PC3, color = Trait)) +
                ggplot2::geom_point(size = size) +
                ggplot2::scale_color_manual(values=custom_cols) +
                ggplot2::labs(x = paste("PC2 (", var_explained["PC2", ], "%)", sep = ""), y = paste("PC3 (", var_explained["PC3", ], "%)", sep = "")) +
                ggplot2::theme_classic() +
                ggplot2::ggtitle("Principal component analysis of samples") +
                ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
        } else {
            stop("Please, specify the PCs to be plotted. One of '1x2', '1x3', or '2x3'.")
        }
    } else {
        if(PCs == "1x2") {
            p <- ggvis::ggvis(pca_df, ~PC1, ~PC2, fill = ~Trait, key := ~SampleID) %>%
                ggvis::layer_points(size := ggvis::input_slider(10, 200, value = 50, label = "Point size")) %>%
                ggvis::add_tooltip(function(df) df$SampleID)
        } else if (PCs == "1x3") {
            p <- ggvis(pca_df, ~PC1, ~PC3, fill = ~Trait, key := ~SampleID) %>%
                ggvis::layer_points(size := ggvis::input_slider(10, 200, value = 50, label = "Point size")) %>%
                ggvis::add_tooltip(function(df) df$SampleID)
        } else if (PCs == "2x3") {
            p <- ggvis(pca_df, ~PC2, ~PC3, fill = ~Trait, key := ~SampleID) %>%
                ggvis::layer_points(size := ggvis::input_slider(10, 200, value = 50, label = "Point size")) %>%
                ggvis::add_tooltip(function(df) df$SampleID)
        } else {
            stop("Please, specify the PCs to be plotted. One of '1x2', '1x3', or '2x3'.")
        }
    }
    p
}


#' Get a vector object with housekeeping genes
#'
#' @param exp Expression dataframe, where rownames are gene IDs and colnames are sample names
#'
#' @return vector object of gene IDs of housekeeping genes
#' @author Fabricio Almeida-Silva
#' @rdname get_HK
#' @export
get_HK <- function(exp) {
    exp=exp
    exp[exp < 1] <- 0 #expression values below 1 are considered as not expressed
    ncols=ncol(exp)
    final.exp <- exp[rowSums(exp > 0) == ncols,]
    final.exp$mean <- rowMeans(final.exp[,1:ncols])
    final.exp$sd <-  apply(final.exp[,1:ncols], 1, sd)
    final.exp$covar <- final.exp$sd/final.exp$mean
    final.exp$max <- apply(final.exp[,1:ncols], 1, max)
    final.exp$min <- apply(final.exp[,1:ncols], 1, min)
    final.exp$MFC <- final.exp$max/final.exp$min
    final.exp$MFC.CoV <- final.exp$MFC * final.exp$covar
    hk <- head(final.exp[order(final.exp$MFC.CoV, decreasing=F),], n = nrow(final.exp)*0.25)
    hk.genes <- rownames(hk)
}

#' Plot expression profile of given genes across samples
#'
#' @param genes Character vector containing a subset of genes from which edges will be extracted. It can be ignored if \code{plot_module} is TRUE.
#' @param exp Data frame containing genes IDs in row names and sample names in column names.
#' @param metadata A 2-column data frame containing sample names in the first column and sample descriptions in the second column.
#' @param plot_module Logical indicating whether to plot a whole module or not. If set to FALSE, \code{genes} must be specified.
#' @param genes_modules Data frame containing genes in column 1 and their corresponding modules in column 2. It is the third element of the output list from \code{exp2net}.
#' @param modulename Character with name of the module to be plotted. To include 2 or more modules, input the names in a character vector.
#'
#' @return A ggplot object showing the expression profile of some genes across all samples.
#' @author Fabricio Almeida-Silva
#' @seealso
#'  \code{\link[reshape2]{melt}}
#' @rdname plot_expression_profile
#' @export
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes_ geom_tile aes geom_line stat_summary theme element_text element_blank element_rect ggtitle
plot_expression_profile <- function(genes, exp, metadata, plot_module = TRUE, genes_modules, modulename) {
    sample_names <- colnames(metadata)[1]
    sample_info <- colnames(metadata)[2]

    if(plot_module == TRUE) {
        genes_in_module <- genes_modules[genes_modules[,2] == modulename, 1]
        filt_exp <- exp[rownames(exp) %in% genes_in_module, ]
        filt_exp$id <- rownames(filt_exp)

        melt_exp <- reshape2::melt(filt_exp, "id",
                                   variable.name = "sample",
                                   value.name = "expression")
    } else {
        melt_exp <- reshape2::melt(filt_exp[genes, ], "id",
                                   variable.name = "sample",
                                   value.name = "expression")
    }
    # order sample info by tissue
    metadata <- metadata[metadata[,1] %in% colnames(filt_exp), ]
    metadata <- metadata[order(metadata[,2]), ]
    metadata[,1] <- factor(metadata[,1], levels = metadata[,1])

    melt_exp[, "sample"] <- factor(melt_exp[, "sample"],
                                   levels = metadata[,1])

    # Background tiles
    background <- mean(melt_exp[, "expression"])

    # Plot expression profiles
    p <- ggplot2::ggplot(melt_exp, ggplot2::aes_(x = ~sample, y = ~expression)) +
        ggplot2::geom_tile(data = metadata, alpha = 0.3, height = Inf,
                           ggplot2::aes(x = get(sample_names), y = background,
                      fill = as.factor(get(sample_info))))

    p <- p + ggplot2::geom_line(ggplot2::aes_(group = ~id), alpha = 0.2,
                       color = "firebrick") +
        ggplot2::stat_summary(ggplot2::aes(group = 1), size = 1, fun = "median", geom = "line")

    p <- p + ggplot2::theme(plot.title=ggplot2::element_text(lineheight=0.8,
                                           face='bold',
                                           colour='black',
                                           size=15),
                   axis.title=ggplot2::element_text(face='bold',
                                           colour='black',
                                           size=15),
                   axis.text.y=ggplot2::element_text(angle=0,
                                            vjust=0.5,
                                            size=8),
                   axis.text.x=ggplot2::element_text(angle=90,
                                            vjust=0.5,
                                            size=6),
                   panel.grid=ggplot2::element_blank(),
                   legend.title=ggplot2::element_blank(),
                   legend.text=ggplot2::element_text(size = 8),
                   legend.background=ggplot2::element_rect(fill='gray90',
                                                  size=0.5,
                                                  linetype='dotted'),
                   legend.position='bottom'
    )

    p <- p + ggplot2::ggtitle(modulename)

    return(p)
}


#' Plot number of genes per module
#'
#' @param genes_and_modules A 2-column data frame with gene IDs in the first column and module names in the second column. Third element of the result list from \code{exp2net}.
#' @return A ggplot object with a barplot of gene number in each module.
#' @seealso
#'  \code{\link[ggpubr]{ggbarplot}}
#' @rdname plot_ngenes_per_module
#' @export
#' @importFrom ggpubr ggbarplot
#' @importFrom ggplot2 theme element_text
plot_ngenes_per_module <- function(genes_and_modules = NULL) {
    frequency_df <- as.data.frame(table(genes_and_modules[,2]),
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


