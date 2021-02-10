
#' Infer gene regulatory network with the Context Likelihood of Relatedness (CLR) algorithm
#'
#' @param exp Expression matrix with gene IDs as row names and samples as column names.
#' @param estimator_clr Entropy estimator to be used. One of "mi.empirical", "mi.mm", "mi.shrink", "mi.sg", "pearson", "spearman", or "kendall". Default: "pearson".
#' @param regulators A character vector of regulators (e.g., transcription factors or miRNAs). All regulators must be included in `exp`.
#' @param remove_zero Logical indicating whether to remove edges whose weight is exactly zero. Default: TRUE
#'
#' @return A gene regulatory network represented as an edge list.
#' @export
#' @rdname grn_clr
#' @importFrom minet build.mim clr
grn_clr <- function(exp, estimator_clr = "pearson",
                    regulators=NULL,
                    remove_zero=TRUE) {

    if(is.null(regulators)) {
        stop("Please, input a character vector of IDs of regulators.")
    }

    # Build Mutual Information matrix and infer GRN
    mi_mat <- minet::build.mim(t(exp), estimator = estimator_clr)
    grn <- minet::clr(mi_mat)

    # Keep only interactions between regulators and targets
    grn <- grn[regulators, !(colnames(grn) %in% regulators)]
    grn_edges <- cormat_to_edgelist(grn)

    # Should we remove edges whose weight is zero?
    if(remove_zero) {
        grn_edges <- grn_edges[grn_edges$Weight != 0, ]
    }

    # Sort weights in decreasing order
    grn_edges <- grn_edges[order(-grn_edges$Weight), ]

    return(grn_edges)
}


#' Infer gene regulatory network with the ARACNE algorithm
#'
#' @param exp Expression matrix with gene IDs as row names and samples as column names.
#' @param estimator_aracne Entropy estimator to be used. One of "mi.empirical", "mi.mm", "mi.shrink", "mi.sg", "pearson", "spearman", or "kendall". Default: "spearman".
#' @param regulators A character vector of regulators (e.g., transcription factors or miRNAs). All regulators must be included in `exp`.
#' @param eps Numeric value indicating the threshold used when removing an edge: for each triplet of nodes (i,j,k), the weakest edge, say (ij), is removed if its weight is below min{(ik),(jk)} - eps. Default: 0.1.
#' @param remove_zero Logical indicating whether to remove edges whose weight is exactly zero. Zero values indicate edges that were removed by ARACNE. Default: TRUE.
#'
#' @return A gene regulatory network represented as an edge list.
#' @export
#' @rdname grn_aracne
#' @importFrom minet aracne
grn_aracne <- function(exp, estimator_aracne = "spearman",
                       regulators=NULL, eps=0.1,
                       remove_zero=TRUE) {

    if(is.null(regulators)) {
        stop("Please, input a character vector of IDs of regulators.")
    }

    # Build Mutual Information matrix and infer GRN
    mi_mat <- minet::build.mim(t(exp), estimator = estimator_aracne)
    grn <- minet::aracne(mi_mat, eps = eps)

    # Keep only interactions between regulators and targets
    grn <- grn[regulators, !(colnames(grn) %in% regulators)]
    grn_edges <- cormat_to_edgelist(grn)

    # Should we remove edges that were removed by ARACNE?
    if(remove_zero) {
        grn_edges <- grn_edges[grn_edges$Weight != 0, ]
    }

    # Sort weights in decreasing order
    grn_edges <- grn_edges[order(-grn_edges$Weight), ]

    return(grn_edges)
}


#' Infer gene regulatory network with GENIE3
#'
#' @param exp Expression matrix with gene IDs as row names and samples as column names.
#' @param regulators A character vector of regulators (e.g., transcription factors or miRNAs). All regulators must be included in `exp`.
#' @param remove_zero Logical indicating whether to remove edges whose weight is exactly zero. Zero values indicate edges that were removed by ARACNE. Default: TRUE.
#' @param ... Additional arguments passed to `GENIE3::GENIE3()`.
#'
#' @return A gene regulatory network represented as an edge list.
#' @importFrom GENIE3 GENIE3
#' @rdname grn_genie3
#' @export
grn_genie3 <- function(exp, regulators = NULL,
                       remove_zero=TRUE, ...) {

    if(is.null(regulators)) {
        stop("Please, input a character vector of IDs of regulators.")
    }

    # Infer GRN
    grn <- GENIE3::GENIE3(as.matrix(exp), regulators = regulators, ...)

    # Keep only interactions between regulators and targets
    grn <- grn[regulators, !(colnames(grn) %in% regulators)]
    grn_edges <- cormat_to_edgelist(grn)

    # Should we remove edges that were removed by ARACNE?
    if(remove_zero) {
        grn_edges <- grn_edges[grn_edges$Weight != 0, ]
    }

    # Sort weights in decreasing order
    grn_edges <- grn_edges[order(-grn_edges$Weight), ]

    return(grn_edges)
}


#' Infer gene regulatory network with multiple algorithms and combine results in a list
#'
#' @param exp Expression matrix with gene IDs as row names and samples as column names.
#' @param estimator_clr Entropy estimator to be used in CLR inference. One of "mi.empirical", "mi.mm", "mi.shrink", "mi.sg", "pearson", "spearman", or "kendall". Default: "pearson".
#' @param estimator_aracne Entropy estimator to be used in ARACNE inference. One of "mi.empirical", "mi.mm", "mi.shrink", "mi.sg", "pearson", "spearman", or "kendall". Default: "spearman".
#' @param regulators A character vector of regulators (e.g., transcription factors or miRNAs). All regulators must be included in `exp`.
#' @param eps Numeric value indicating the threshold used when removing an edge: for each triplet of nodes (i,j,k), the weakest edge, say (ij), is removed if its weight is below min{(ik),(jk)} - eps. Default: 0.
#' @param remove_zero Logical indicating whether to remove edges whose weight is exactly zero. Zero values indicate edges that were removed by ARACNE. Default: TRUE.
#' @param ... Additional arguments passed to `GENIE3::GENIE3()`.
#'
#' @return A list of data frames representing edge lists. Each list element is an edge list for a specific method.
#' @rdname grn_combined
#' @export
grn_combined <- function(exp, regulators = NULL,
                         eps=0,
                         estimator_aracne = "spearman",
                         estimator_clr = "pearson",
                         remove_zero=TRUE, ...) {

    genie3 <- grn_genie3(exp, regulators, remove_zero=remove_zero, ...)

    aracne <- grn_aracne(exp, regulators = regulators, eps=eps,
                         estimator_aracne = estimator_aracne, remove_zero=remove_zero)

    clr <- grn_clr(exp, regulators = regulators,
                   estimator_clr = estimator_clr, remove_zero=remove_zero)

    res_list <- list(genie3=genie3, aracne=aracne, clr=clr)
}

#' Rank edge weights for GRNs and calculate average across different methods
#'
#' @param list_edges List containing edge lists as returned by the function \code{grn_combined}.
#'
#' @return Edge list containing regulator, target and mean rank from all algorithms.
#' @rdname grn_average_rank
#' @export
grn_average_rank <- function(list_edges) {

    # Add ranks to each edge list
    ranked_list <- lapply(list_edges, function(x) {
        x$Rank <- seq_len(nrow(x))
        x$Weight <- NULL
        return(x)
    })

    # Reduce list to a single data frame
    df_ranks <- Reduce(function(x, y) merge(x, y, by=c("Node1", "Node2")), ranked_list)
    colnames(df_ranks) <- c("Regulator", "Target", "Rank1", "Rank2", "Rank3")

    # Calculate mean and order from the highest mean rank to the lowest
    df_ranks$mean_rank <- apply(df_ranks[, 3:ncol(df_ranks)], 1, mean)
    df_ranks <- df_ranks[, c("Regulator", "Target", "mean_rank")]
    df_ranks <- df_ranks[order(df_ranks$mean_rank), ]
    return(df_ranks)
}



#' Filter a gene regulatory network based on optimal scale-free topology fit
#'
#' @param edgelist A gene regulatory network represented as an edge list.
#' @param nsplit Number of groups in which the edge list will be split. Default: 10.
#'
#' @details The edge list will be split in n groups and the scale-free topology fit will be tested for each subset of the edge list.
#' For instance, if an edge list of 10000 rows is used as input, the function will test SFT fit for the top 1000 edges, then top 2000 edges, an so on up to the whole edge list.
#' @return The edge list that best fits the scale-free topology.
#' @export
#' @importFrom igraph graph_from_data_frame degree
#' @importFrom BiocParallel bplapply
#' @importFrom WGCNA scaleFreeFitIndex
#' @importFrom ggpubr ggline
#' @importFrom ggplot2 theme element_text
grn_filter <- function(edgelist, nsplit=10) {

    # Split edge list into n data frames and calculate degree
    n_edges <- seq_len(nrow(edgelist))
    cutpoints <- round(unname(quantile(n_edges, probs = seq(0, 1, 1/nsplit))))[-1]

    list_degree <- BiocParallel::bplapply(cutpoints, function(x) {
        filt_edges <- edgelist[1:x, 1:2]
        graph <- igraph::graph_from_data_frame(filt_edges, directed=TRUE)
        degree <- igraph::degree(graph, mode = "out")
    })

    # Calculate scale-free topology fit for the degree list
    sft.rsquared <- unlist(BiocParallel::bplapply(list_degree, function(x) {
        return(WGCNA::scaleFreeFitIndex(x)$Rsquared.SFT)
        }))
    max.index <- which.max(sft.rsquared)

    # Plot scale-free topology fit for r values
    plot.data <- data.frame(x=cutpoints, y=sft.rsquared, stringsAsFactors = FALSE)
    plot <- ggpubr::ggline(plot.data, x = "x", y = "y", size=2,
                           color="firebrick",
                           xlab = "Number of top edges considered",
                           xtickslab.rt = 45,
                           ylab = expression(paste("Scale-free topology fit - ", R^{2})),
                           title = "Scale-free topology fit for top n edges", font.title = c(13, "bold")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    print(plot)

    optimal_cutoff <- cutpoints[max.index]
    message("The top number of edges that best fits the scale-free topology is ", optimal_cutoff)

    edgelist <- edgelist[1:optimal_cutoff, 1:2]
    return(edgelist)
}


#' Get hubs for gene regulatory network
#'
#' @param edgelist A gene regulatory network represented as an edge list.
#' @param top_percentile Numeric from 0 to 1 indicating the percentage of genes with the highest degree to consider hubs. Default: 0.1.
#' @param top_n Numeric indicating the number of genes with the highest degree to consider hubs.
#'
#' @return A data frame with gene ID in the first column and out degree in the second column.
#' @export
#' @rdname get_hubs_grn
#' @importFrom igraph graph_from_data_frame degree
get_hubs_grn <- function(edgelist, top_percentile = 0.1, top_n = NULL) {

    # Calculate degree
    graph <- igraph::graph_from_data_frame(edgelist, directed = TRUE)
    degree <- sort(igraph::degree(graph, mode="out"), decreasing = TRUE)

    # Find hubs
    degree_df <- data.frame(row.names=1:length(degree),
                            Gene=names(degree),
                            Degree=degree, stringsAsFactors = FALSE)

    if(is.null(top_n)) {
        nrows <- nrow(degree_df) * top_percentile
        hubs <- degree_df[1:nrows, ]
    } else {
        hubs <- degree_df[1:top_n, ]
    }

    return(hubs)
}













