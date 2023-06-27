
#' Helper to handle SummarizedExperiment or expression data frame as input
#'
#' @param exp Expression data as a data frame or a SummarizedExperiment object
#'
#' @return If exp is a SummarizedExperiment object, it will return `assay(exp)`.
#' Otherwhise, it will simply return exp as it is.
#' @noRd
#' @importFrom SummarizedExperiment assay
handleSE <- function(exp) {
    if(is(exp, "SummarizedExperiment")) {
        fexp <- SummarizedExperiment::assay(exp)
    } else {
        fexp <- exp
    }
    return(fexp)
}

#' Replace content of a SummarizedExperiment object based on filtered expression data frame
#'
#' @param exp Expression data frame with genes IDs in row names and samples
#' in column names.
#' @param SE Original SummarizedExperiment object.
#' @return A SummarizedExperiment object
#' @noRd
#' @importFrom SummarizedExperiment SummarizedExperiment colData assay
exp2SE <- function(exp, SE) {

    # Get overlaps between samples in both sets
    SE_cols <- rownames(SummarizedExperiment::colData(SE))
    filt_exp_cols <- colnames(exp)
    overlap <- which(SE_cols %in% filt_exp_cols)

    # Modify original SE based on filtered expression data frame
    SE_final <- SummarizedExperiment::SummarizedExperiment(
        assays = exp,
        colData = SummarizedExperiment::colData(SE)[overlap, , drop = FALSE]
    )

    return(SE_final)
}


#' Generate custom color palette
#'
#' @param pal Numeric specifying palette number, from 1 to 3.
#'
#' @return Character vector of custom color palette with 20 colors
#' @noRd
custom_palette <- function(pal = 1) {
    pal1 <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF",
              "#9467BDFF", "#8C564BFF", "#E377C2FF", "#7F7F7FFF",
              "#BCBD22FF", "#17BECFFF", "#AEC7E8FF", "#FFBB78FF",
              "#98DF8AFF", "#FF9896FF", "#C5B0D5FF", "#C49C94FF",
              "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF")

    pal2 <- c("#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF",
              "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF",
              "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF",
              "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", "#C6DBEFFF",
              "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF")

    pal3 <- c("#393B79FF", "#637939FF", "#8C6D31FF", "#843C39FF",
              "#7B4173FF", "#5254A3FF", "#8CA252FF", "#BD9E39FF",
              "#AD494AFF", "#A55194FF", "#6B6ECFFF", "#B5CF6BFF",
              "#E7BA52FF", "#D6616BFF", "#CE6DBDFF", "#9C9EDEFF",
              "#CEDB9CFF", "#E7CB94FF", "#E7969CFF", "#DE9ED6FF")

    l <- list(pal1, pal2, pal3)
    l_final <- l[[pal]]
    return(l_final)
}

#' Map levels of metadata variables to colors for plotting
#'
#' @param col_metadata A data frame with column or row metadata. If column
#' metadata is passed, row names must contain sample names. If row metadata is
#' passed, row names must contain gene names.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{metadata}{A metadata data frame as in \strong{col_metadata}, but
#'                   with rows sorted by levels of every column.}
#'   \item{colors}{A list of named character vectors containing the mapping
#'                 between levels of metadata variables and colors.}
#' }
#'
#' @noRd
#' @importFrom stats setNames
metadata2colors <- function(col_metadata) {

    coldata <- NA
    colors <- NA

    if(is.data.frame(col_metadata)) {
        # Get variable names in metadata
        col_names <- names(col_metadata)
        if(length(col_names) > 3) {
            stop("Maximum number of columns for row and sample metadata is 3.")
        }

        # Sort elements in all columns
        coldata <- col_metadata[do.call(order, col_metadata), , drop = FALSE]

        # Create a list of named vectors with variable levels and colors
        colors <- lapply(seq_along(col_names), function(x) {
            levels <- unique(coldata[, x])
            cols <- setNames(custom_palette(x)[seq_along(levels)], levels)
            return(cols)
        })
        names(colors) <- col_names
    }

    # Return results as a list
    results <- list(
        metadata = coldata,
        colors = colors
    )

    return(results)
}

#' Construct parameters to plot heatmap
#'
#' @param exp A gene expression matrix with gene IDs in row names and
#' samples in column names.
#' @param palette Character indicating an RColorBrewer palette to use.
#' @param heatmap_type Character indicating which heatmap type to plot.
#' One of "samplecor" (for pairwise sample correlations) or "expr" (for
#' gene expression).
#' @param cor_method Character indicating which correlation method to use
#' in the \code{cor()} function. Only meaningful if parameter
#' \strong{heatmap_type} equals "samplecor".
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{pal}{Character, name of the RColorBrewer palette to use.}
#'   \item{mat}{Matrix to use to construct the heatmap.}
#'   \item{title}{Character, title of the heatmap.}
#'   \item{name}{Character, matrix name. This is what goes to the title of
#'               matrix's main legend.}
#' }
#'
#' @noRd
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats cor
#'
heatmap_attributes <- function(exp, palette = NULL, heatmap_type = "samplecor",
                               cor_method = "pearson") {

    # Create objects with attributes `pal`, `mat`, `title`, and `name`
    if(heatmap_type == "samplecor") {
        pal <- "Blues"
        mat <- cor(exp, method = cor_method)
        title <- "Pairwise correlations between samples"
        name <- "Correlation"
    } else if(heatmap_type == "expr") {
        pal <- "YlOrRd"
        mat <- exp
        title <- "Gene expression"
        name <- "Expression"
    } else {
        stop("Invalid argument to `type`. Pick one of 'samplecor' or 'expr'.")
    }
    if(!is.null(palette)) { pal <- palette }
    pal <- colorRampPalette(RColorBrewer::brewer.pal(9, pal))(100)

    # Return results as a list
    results <- list(
        pal = pal,
        mat = mat,
        title = title,
        name = name
    )
    return(results)
}

#' Extract row and column metadata from `SummarizedExperiment` objects
#'
#' @param se A `SummarizedExperiment` object.
#' @param rowdata_cols Columns to use from the rowData element of the
#' `SummarizedExperiment` object. It can be either a character vector
#' of column names or a numeric vector of column indices.
#' By default, all columns are used.
#' @param coldata_cols Columns to use from the colData element of the
#' `SummarizedExperiment` object. It can be either a character vector
#' of column names or a numeric vector of column indices.
#' By default, all columns are used.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{rowdata}{A data frame of row metadata containing only the selected
#'                  columns.}
#'   \item{coldata}{A data frame of column metadata containing only the
#'                  selected columns.}
#' }
#'
#' @noRd
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment colData
se2metadata <- function(se, rowdata_cols = NULL, coldata_cols = NULL) {

    final_rowdata <- NA
    final_coldata <- NA

    # Extract row and column metadata
    rowdata <- SummarizedExperiment::rowData(se)
    coldata <- SummarizedExperiment::colData(se)

    # Extract user-defined columns from metadata
    ## Row metadata
    if(ncol(rowdata) > 0) {
        r_cols <- rowdata_cols
        if(is.null(rowdata_cols)) { r_cols <- seq_along(rowdata) }

        final_rowdata <- as.data.frame(rowdata[, r_cols, drop = FALSE])
    }

    if(ncol(coldata) > 0) {
        c_cols <- coldata_cols
        if(is.null(coldata_cols)) { c_cols <- seq_along(coldata) }
        final_coldata <- as.data.frame(coldata[, c_cols, drop = FALSE])
    }

    # Return resuls as a list
    metadata_list <- list(
        rowdata = final_rowdata,
        coldata = final_coldata
    )

    return(metadata_list)
}


#' Get model matrix for metadata variables
#'
#' @param metadata A data frame of column metadata with sample names
#' in row names.
#' @param column_idx Column to use to create the model matrix.
#'
#' @details If the variable is numeric (continuous or discrete),
#' the model matrix is created by simply subsetting the column
#' indicated in \strong{column_idx}. If the variable if categorical,
#' a dummy model matrix is created (without an intercept).
#'
#' @return A data frame with the model matrix to use
#' in \code{module_trait_cor()}.
#'
#' @noRd
#' @importFrom stats model.matrix as.formula
get_model_matrix <- function(metadata, column_idx = 1) {

    # By default, assuming variable is numeric (continuous or discrete)
    mat <- metadata[, column_idx, drop = FALSE]

    if(!is.numeric(metadata[, column_idx])) {

        # Create formula (format: `~var + 0`)
        var_name <- names(metadata)[column_idx]
        formula <- as.formula(paste("~", var_name, "+ 0"))

        # Get model matrix
        mat <- as.data.frame(model.matrix(formula, metadata))
        names(mat) <- gsub(var_name, "", names(mat))
    }

    return(mat)
}


#' Wrapper to calculate correlation matrix and adjacency matrices
#'
#' @param cor_method Correlation method.
#' @param norm.exp Gene expression data frame.
#' @param SFTpower SFT power calculated with \code{SFT_fit}.
#' @param net_type Network type. One of 'signed', 'signed hybrid' or 'unsigned'.
#'
#' @return List containing the correlation matrix and the adjacency matrix.
#' @noRd
#' @importFrom WGCNA adjacency.fromSimilarity bicor
calculate_cor_adj <- function(cor_method, norm.exp, SFTpower,
                              net_type) {

    if(cor_method == "pearson") {
        cor_matrix <- cor(t(norm.exp), method = "pearson")
        adj_matrix <- WGCNA::adjacency.fromSimilarity(
            cor_matrix, power = SFTpower, type = net_type
        )
    } else if(cor_method == "spearman") {
        cor_matrix <- cor(t(norm.exp), use="p", method = "spearman")
        adj_matrix <- WGCNA::adjacency.fromSimilarity(
            cor_matrix, power = SFTpower, type = net_type
        )
    } else if (cor_method == "biweight") {
        cor_matrix <- WGCNA::bicor(t(norm.exp), maxPOutliers = 0.1)
        adj_matrix <- WGCNA::adjacency.fromSimilarity(
            cor_matrix, power = SFTpower, type = net_type
        )
    } else {
        stop("Please, specify a correlation method. One of 'spearman', 'pearson' or 'biweight'.")
    }
    results <- list(cor = cor_matrix, adj = adj_matrix)
    return(results)
}


#' Wrapper to assign TOM type
#'
#' @param net_type Network type. One of 'signed', 'signed hybrid' or 'unsigned'.
#'
#' @return Character of TOM type
#' @noRd
get_TOMtype <- function(net_type) {
    if(net_type == "signed hybrid") {
        TOMType <- "signed"
    } else if(net_type == "signed") {
        TOMType <- "signed Nowick"
    } else {
        TOMType <- "unsigned"
    }
    return(TOMType)
}


#' Transform a correlation matrix to an edge list
#'
#' @param matrix Symmetrical correlation matrix.
#'
#' @return A 2-column data frame containing node 1, node 2 and edge weight.
#' @export
#' @rdname cormat_to_edgelist
#' @examples
#' data(filt.se)
#' cor_mat <- cor(t(SummarizedExperiment::assay(filt.se)))
#' edgelist <- cormat_to_edgelist(cor_mat)
cormat_to_edgelist <- function(matrix) {

    edgelist <- matrix
    edgelist[lower.tri(edgelist, diag=TRUE)] <- NA
    edgelist <- na.omit(data.frame(as.table(edgelist), stringsAsFactors=FALSE))
    colnames(edgelist) <- c("Node1", "Node2", "Weight")
    edgelist$Node1 <- as.character(edgelist$Node1)
    edgelist$Node2 <- as.character(edgelist$Node2)
    edgelist$Weight <- as.numeric(edgelist$Weight)
    return(edgelist)
}


#' Check scale-free topology fit for a given network
#'
#' @param edgelist Edge list as a data frame containing node 1,
#' node 2 and edge weight.
#' @param net_type Type of biological network. One of "gcn", "grn", or "ppi".
#' Default: gcn.
#'
#' @return A list with SFT fit statistics and a message indicating if
#' the network is scale-free.
#' @rdname check_SFT
#' @export
#' @importFrom igraph graph_from_data_frame as_adjacency_matrix fit_power_law
#' @examples
#' set.seed(1)
#' exp <- t(matrix(rnorm(10000), ncol=1000, nrow=200))
#' rownames(exp) <- paste0("Gene", 1:nrow(exp))
#' colnames(exp) <- paste0("Sample", 1:ncol(exp))
#' cormat <- cor(t(exp))
#' edges <- cormat_to_edgelist(cormat)
#' edges <- edges[abs(edges$Weight) > 0.10, ]
#' check_SFT(edges)
check_SFT <- function(edgelist, net_type = "gcn") {

    # Calculate degree of the resulting graph
    if(net_type == "gcn") {
        graph <- igraph::graph_from_data_frame(edgelist, directed = FALSE)
        adj <- igraph::as_adjacency_matrix(graph, sparse = FALSE)
        diag(adj) <- 0
        degree <- apply(adj, 1, sum, na.rm = TRUE)
    } else if(net_type == "grn") {
        graph <- igraph::graph_from_data_frame(edgelist, directed = TRUE)
        degree <- igraph::degree(graph, mode = "out")
    } else if(net_type == "ppi") {
        graph <- igraph::graph_from_data_frame(edgelist, directed = FALSE)
        degree <- igraph::degree(graph)
    } else {
        stop("Invalid network type. Please, input one of 'gcn', 'grn', or 'ppi'.")
    }

    # Test for scale-free topology fit
    test <- igraph::fit_power_law(degree)
    if(test$KS.p < 0.05) {
        message("At the 95% confidence level for the Kolmogorov-Smirnov statistic, your graph does not fit the scale-free topology. P-value:", test$KS.p)
    } else {
        message("Your graph fits the scale-free topology. P-value:", test$KS.p)
    }

    return(test)
}


#' Helper to handle list of SummarizedExperiment objects as input
#'
#' @param exp List of data frames or SummarizedExperiment objects.
#'
#' @return If exp is a list of SummarizedExperiment objects,
#' it will return a list of data frames. Otherwise, it will simply return exp as it is.
#' @noRd
#' @importFrom SummarizedExperiment assay
handleSElist <- function(exp) {

    if(is(exp[[1]], "SummarizedExperiment")) {
        list <- lapply(exp, handleSE)
    } else {
        list <- exp
    }
    return(list)
}


#' Helper function to handle metadata input for consensus modules identification
#'
#' @param exp List of data frames or SummarizedExperiment objects.
#' @param metadata A data frame containing sample names in row names and sample annotation in the first column.
#' @return Data frame of metadata for all expression sets.
#' @noRd
#' @importFrom SummarizedExperiment colData
handle_metadata <- function(exp, metadata) {

    if(is(exp[[1]], "SummarizedExperiment")) {
        metadata <- Reduce(rbind, lapply(exp, function(x) {
            meta <- as.data.frame(SummarizedExperiment::colData(x))
            return(meta)
        }))
    } else {
        metadata <- metadata
    }
    return(metadata)
}



#' Convert p-values in matrix to symbols
#'
#' @param matrix Matrix of p-values.
#'
#' @return Matrix of symbols.
#' @noRd
pval2symbol <- function(matrix) {
    modtraitsymbol <- matrix
    modtraitsymbol[modtraitsymbol < 0.001] <- "***"
    modtraitsymbol[modtraitsymbol >= 0.001 & modtraitsymbol < 0.01] <- "**"
    modtraitsymbol[modtraitsymbol >= 0.01 & modtraitsymbol < 0.05] <- "*"
    modtraitsymbol[modtraitsymbol >= 0.05] <- ""
    return(modtraitsymbol)
}















