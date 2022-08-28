

#' Infer gene regulatory network with one of three algorithms
#'
#' The available algorithms are Context Likelihood of Relatedness (CLR),
#' ARACNE, or GENIE3.
#'
#' @param exp A gene expression data frame with genes in row names and
#' samples in column names or a `SummarizedExperiment` object.
#' @param regulators A character vector of regulators
#' (e.g., transcription factors or miRNAs). All regulators must be
#' included in `exp`.
#' @param method GRN inference algorithm to be used. One of "clr", "aracne",
#' or "genie3".
#' @param estimator_clr Entropy estimator to be used. One of "mi.empirical",
#' "mi.mm", "mi.shrink", "mi.sg", "pearson", "spearman", or "kendall".
#' Default: "pearson".
#' @param estimator_aracne Entropy estimator to be used. One of "mi.empirical",
#' "mi.mm", "mi.shrink", "mi.sg", "pearson", "spearman", or "kendall".
#' Default: "spearman".
#' @param eps Numeric value indicating the threshold used when removing
#' an edge: for each triplet of nodes (i,j,k), the weakest edge, say (ij),
#' is removed if its weight is below min{(ik),(jk)} - eps. Default: 0.1.
#' @param remove_zero Logical indicating whether to remove edges whose weight
#' is exactly zero. Default: TRUE
#' @param ... Additional arguments passed to `GENIE3::GENIE3()`.
#' @return A gene regulatory network represented as an edge list.
#' @rdname grn_infer
#' @export
#' @importFrom minet build.mim clr aracne
#' @importFrom GENIE3 GENIE3
#' @examples
#' data(filt.se)
#' tfs <- sample(rownames(filt.se), size=20, replace=FALSE)
#' clr <- grn_infer(filt.se, method = "clr", regulators=tfs)
#' aracne <- grn_infer(filt.se, method = "aracne", regulators=tfs)
#' # only 2 trees for demonstration purposes
#' genie3 <- grn_infer(filt.se, method = "genie3", regulators=tfs, nTrees=2)
grn_infer <- function(exp, regulators=NULL,
                      method = c("clr", "aracne", "genie3"),
                      estimator_clr = "pearson",
                      estimator_aracne = "spearman",
                      eps = 0.1,
                      remove_zero = TRUE, ...) {
    exp <- handleSE(exp)
    if(is.null(regulators)) {
        stop("Please, input a character vector of IDs of regulators.")
    }
    regulators <- regulators[regulators %in% rownames(exp)]

    if(method == "clr") {
        mi_mat <- minet::build.mim(t(exp), estimator = estimator_clr)
        grn <- minet::clr(mi_mat)
    } else if(method == "aracne") {
        # Build Mutual Information matrix and infer GRN
        mi_mat <- minet::build.mim(t(exp), estimator = estimator_aracne)
        grn <- minet::aracne(mi_mat, eps = eps)
    } else if(method == "genie3") {
        grn <- GENIE3::GENIE3(as.matrix(exp), regulators = regulators, ...)
    } else {
        stop("Method must be one of 'clr', 'aracne', or 'genie3'.")
    }

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
#' @param exp A gene expression data frame with genes in row names and
#' samples in column names or a `SummarizedExperiment` object.
#' @param regulators A character vector of regulators
#' (e.g., transcription factors or miRNAs). All regulators must be
#' included in `exp`.
#' @param eps Numeric value indicating the threshold used when removing
#' an edge: for each triplet of nodes (i,j,k), the weakest edge, say (ij),
#' is removed if its weight is below min{(ik),(jk)} - eps. Default: 0.1.
#' @param estimator_clr Entropy estimator to be used in CLR inference.
#' One of "mi.empirical", "mi.mm", "mi.shrink", "mi.sg", "pearson", "spearman",
#' or "kendall". Default: "pearson".
#' @param estimator_aracne Entropy estimator to be used in ARACNE inference.
#' One of "mi.empirical", "mi.mm", "mi.shrink", "mi.sg", "pearson", "spearman",
#' or "kendall". Default: "spearman".
#' @param remove_zero Logical indicating whether to remove edges
#' whose weight is exactly zero. Zero values indicate edges that were
#' removed by ARACNE. Default: TRUE.
#' @param ... Additional arguments passed to `GENIE3::GENIE3()`.
#' @return A list of data frames representing edge lists. Each list element
#' is an edge list for a specific method.
#' @rdname grn_combined
#' @export
#' @examples
#' data(filt.se)
#' tfs <- sample(rownames(filt.se), size=50, replace=FALSE)
#' grn_list <- grn_combined(filt.se, regulators=tfs, nTrees=2)
grn_combined <- function(exp, regulators = NULL,
                         eps = 0.1,
                         estimator_aracne = "spearman",
                         estimator_clr = "pearson",
                         remove_zero = TRUE, ...) {
    regulators <- regulators[regulators %in% rownames(exp)]
    methods <- c("genie3", "aracne", "clr")
    res_list <- lapply(methods, function(x) {
        grns <- grn_infer(exp, regulators = regulators, method = x,
                          estimator_aracne = estimator_aracne,
                          estimator_clr = estimator_clr,
                          remove_zero = remove_zero, ...)
        return(grns)
    })
    names(res_list) <- methods
    return(res_list)
}

#' Rank edge weights for GRNs and calculate average across different methods
#'
#' @param list_edges List containing edge lists as returned by
#' the function \code{grn_combined}.
#'
#' @return Edge list containing regulator, target and mean rank from
#' all algorithms.
#' @rdname grn_average_rank
#' @export
#' @examples
#' data(filt.se)
#' tfs <- sample(rownames(filt.se), size=50, replace=FALSE)
#' grn_list <- grn_combined(filt.se, regulators=tfs, nTrees=2)
#' ranked_grn <- grn_average_rank(grn_list)
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
#' @param nsplit Number of groups in which the edge list will be split.
#' Default: 10.
#' @param bp_param BiocParallel back-end to be used.
#' Default: BiocParallel::SerialParam()
#'
#' @details The edge list will be split in n groups and the scale-free
#' topology fit will be tested for each subset of the edge list.
#' For instance, if an edge list of 10000 rows is used as input,
#' the function will test SFT fit for the top 1000 edges, then top 2000 edges,
#' and so on up to the whole edge list.
#' @return The edge list that best fits the scale-free topology.
#' @export
#' @rdname grn_filter
#' @importFrom igraph graph_from_data_frame degree
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom WGCNA scaleFreeFitIndex
#' @importFrom ggplot2 theme element_text geom_point geom_line theme_bw
#' @examples
#' data(filt.se)
#' tfs <- sample(rownames(filt.se), size=50, replace=FALSE)
#' grn_list <- grn_combined(filt.se, regulators=tfs, nTrees=2)
#' ranked_grn <- grn_average_rank(grn_list)
#' # split in only 2 groups for demonstration purposes
#' filtered_edges <- grn_filter(ranked_grn, nsplit=2)
grn_filter <- function(edgelist, nsplit = 10,
                       bp_param = BiocParallel::SerialParam()) {

    # Split edge list into n data frames and calculate degree
    n_edges <- seq_len(nrow(edgelist))
    cutpoints <- round(unname(quantile(n_edges, probs = seq(0, 1, 1/nsplit))))[-1]

    list_degree <- BiocParallel::bplapply(cutpoints, function(x) {
        filt_edges <- edgelist[seq_len(x), c(1,2)]
        graph <- igraph::graph_from_data_frame(filt_edges, directed=TRUE)
        degree <- igraph::degree(graph, mode = "out")
    }, BPPARAM = bp_param)

    # Calculate scale-free topology fit for the degree list
    sft.rsquared <- unlist(BiocParallel::bplapply(list_degree, function(x) {
        return(WGCNA::scaleFreeFitIndex(x)$Rsquared.SFT)
    }, BPPARAM = bp_param))
    max.index <- which.max(sft.rsquared)

    # Plot scale-free topology fit for r values
    plot.data <- data.frame(x = cutpoints, y = sft.rsquared)

    plot <- ggplot(plot.data, aes_(x = ~x, y = ~y, group = 1)) +
        geom_point(color = "firebrick", size = 4) +
        geom_line(color = "firebrick", size = 2) +
        labs(
            x = "Number of top edges considered",
            y = expression(paste("Scale-free topology fit - ", R^{2})),
            title = "Scale-free topology fit for given r values"
        ) +
        theme_bw()
    print(plot)

    optimal_cutoff <- cutpoints[max.index]
    message("The top number of edges that best fits the scale-free topology is ", optimal_cutoff)

    edgelist <- edgelist[seq_len(optimal_cutoff), c(1,2)]
    return(edgelist)
}


#' Infer gene regulatory network from expression data
#'
#' @param exp A gene expression data frame with genes in row names and
#' samples in column names or a `SummarizedExperiment` object.
#' @param estimator_clr Entropy estimator to be used in CLR inference.
#' One of "mi.empirical", "mi.mm", "mi.shrink", "mi.sg", "pearson",
#' "spearman", or "kendall". Default: "pearson".
#' @param estimator_aracne Entropy estimator to be used in ARACNE inference.
#' One of "mi.empirical", "mi.mm", "mi.shrink", "mi.sg", "pearson", "spearman",
#' or "kendall". Default: "spearman".
#' @param regulators A character vector of regulators
#' (e.g., transcription factors or miRNAs). All regulators must be
#' included in `exp`.
#' @param eps Numeric value indicating the threshold used when
#' removing an edge: for each triplet of nodes (i,j,k), the weakest edge,
#' say (ij), is removed if its weight is below min{(ik),(jk)} - eps. Default: 0.
#' @param remove_zero Logical indicating whether to remove edges whose
#' weight is exactly zero. Zero values indicate edges that were
#' removed by ARACNE. Default: TRUE.
#' @param nsplit Number of groups in which the edge list will be split.
#' Default: 10.
#' @param ... Additional arguments passed to `GENIE3::GENIE3()`.
#' @details
#' This function infers GRNs with ARACNE, GENIE3 and CLR, ranks correlation
#' weights for each GRN and calculates the average rank for each edge.
#' Then, the resulting GRN is filtered to keep the top n edges that lead
#' to the optimal scale-free topology fit.
#'
#' @rdname exp2grn
#' @export
#' @return A filtered edge list with regulators in the first column and
#' targets in the second column.
#' @examples
#' data(filt.se)
#' tfs <- sample(rownames(filt.se), size=50, replace=FALSE)
#' # Test with small number of trees for demonstration purpose
#' grn <- exp2grn(filt.se, regulators = tfs, nTrees=2, nsplit=2)
exp2grn <- function(exp, regulators = NULL,
                    eps=0,
                    estimator_aracne = "spearman",
                    estimator_clr = "pearson",
                    remove_zero=TRUE,
                    nsplit=10,
                    ...) {
    grn_list <- grn_combined(exp, regulators = regulators,
                             eps = eps,
                             estimator_aracne = estimator_aracne,
                             estimator_clr = estimator_clr,
                             remove_zero = remove_zero, ...)

    grn_ranks <- grn_average_rank(grn_list)
    final_grn <- grn_filter(grn_ranks, nsplit=nsplit)
    return(final_grn)
}

#' Get hubs for gene regulatory network
#'
#' @param edgelist A gene regulatory network represented as an edge list.
#' @param top_percentile Numeric from 0 to 1 indicating the percentage of genes
#' with the highest degree to consider hubs. Default: 0.1.
#' @param top_n Numeric indicating the number of genes with the highest degree
#' to consider hubs.
#' @param return_degree Logical indicating whether to return a data frame of
#' degree for all genes. If TRUE, the function will return a list instead of
#' a data frame. Default: FALSE.
#' @param ranked Logical indicating whether to treat third column of the
#' edge list (edge weights) as ranked values. Ignored if the edge list only
#' contains 2 columns. Default: TRUE.
#' @return A data frame with gene ID in the first column and out degree in
#' the second column or a list of two data frames with hubs and degree for
#' all genes, respectively.
#' @export
#' @rdname get_hubs_grn
#' @importFrom igraph graph_from_data_frame degree
#' @examples
#' data(filt.se)
#' tfs <- sample(rownames(filt.se), size=50, replace=FALSE)
#' grn_list <- grn_combined(filt.se, regulators=tfs, nTrees=2)
#' ranked_grn <- grn_average_rank(grn_list)
#' # split in only 2 groups for demonstration purposes
#' filtered_edges <- grn_filter(ranked_grn, nsplit=2)
#' hubs <- get_hubs_grn(filtered_edges)
get_hubs_grn <- function(edgelist, top_percentile = 0.1, top_n = NULL,
                         return_degree=FALSE, ranked=TRUE) {

    # Calculate degree
    if(ncol(edgelist) != 2) {
        if(ranked) {
            edgelist[,3] <- 10 / edgelist[,3]
        }
    }
    graph <- igraph::graph_from_data_frame(edgelist, directed = TRUE)
    degree <- sort(igraph::degree(graph, mode="out"), decreasing = TRUE)

    # Find hubs
    degree_df <- data.frame(row.names=seq_along(degree),
                            Gene=names(degree),
                            Degree=degree, stringsAsFactors = FALSE)

    if(is.null(top_n)) {
        nrows <- nrow(degree_df) * top_percentile
        hubs <- degree_df[seq_len(nrows), ]
    } else {
        hubs <- degree_df[seq_len(top_n), ]
    }
    results <- hubs
    if(return_degree) {
        results <- list(Hubs=hubs,
                        Degree=degree_df)
    }
    return(results)
}


#' Get hubs for protein-protein interaction network
#'
#' @param edgelist A protein-protein interaction network represented
#' as an edge list.
#' @param top_percentile Numeric from 0 to 1 indicating the percentage of
#' proteins with the highest degree to consider hubs. Default: 0.1.
#' @param top_n Numeric indicating the number of proteins with the highest
#' degree to consider hubs.
#' @param return_degree Logical indicating whether to return a data frame
#' of degree for all proteins. If TRUE, the function will return a list
#' instead of a data frame. Default: FALSE.
#' @return A data frame with protein ID in the first column and degree
#' in the second column or a list of two data frames with hubs and degree
#' for all genes, respectively.
#' @export
#' @rdname get_hubs_grn
#' @importFrom igraph graph_from_data_frame degree
#' @examples
#' ppi_edges <- igraph::get.edgelist(igraph::barabasi.game(n=500, directed=FALSE))
#' hubs <- get_hubs_ppi(ppi_edges, return_degree = TRUE)
get_hubs_ppi <- function(edgelist, top_percentile = 0.1, top_n = NULL,
                         return_degree=FALSE) {

    # Calculate degree
    graph <- igraph::graph_from_data_frame(edgelist, directed = FALSE)
    degree <- sort(igraph::degree(graph), decreasing = TRUE)

    # Find hubs
    degree_df <- data.frame(row.names=seq_along(degree),
                            Protein=names(degree),
                            Degree=degree, stringsAsFactors = FALSE)

    if(is.null(top_n)) {
        nrows <- nrow(degree_df) * top_percentile
        hubs <- degree_df[seq_len(nrows), ]
    } else {
        hubs <- degree_df[seq_len(top_n), ]
    }
    results <- hubs
    if(return_degree) {
        results <- list(Hubs=hubs,
                        Degree=degree_df)
    }
    return(results)
}












