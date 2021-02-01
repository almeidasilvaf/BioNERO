#' Transform a correlation matrix to an edge list
#'
#' @param matrix Symmetrical correlation matrix.
#'
#' @return A 2-column data frame containing node 1, node 2 and edge weight
#' @export
#' @rdname cormat_to_edgelist
#' @export
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


#' Infer gene regulatory network with the Context Likelihood of Relatedness (CLR) algorithm
#'
#' @param exp Expression matrix with gene IDs as row names and samples as column names.
#' @param estimator Entropy estimator to be used. One of "mi.empirical", "mi.mm", "mi.shrink", "mi.sg", "pearson", "spearman", or "kendall". Default: "pearson".
#'
#' @return A gene regulatory network represented as an edge list.
#' @export
#' @rdname grn_clr
#' @importFrom minet build.mim clr
grn_clr <- function(exp, estimator = "pearson") {
    mi_mat <- minet::build.mim(t(exp), estimator = estimator)
    grn <- minet::clr(mi_mat)
    grn_edges <- cormat_to_edgelist(grn)
    grn_edges <- grn_edges[order(-grn_edges[,3]), ]
    return(grn_edges)
}


#' Infer gene regulatory network with the ARACNE algorithm
#'
#' @param exp Expression matrix with gene IDs as row names and samples as column names.
#' @param estimator Entropy estimator to be used. One of "mi.empirical", "mi.mm", "mi.shrink", "mi.sg", "pearson", "spearman", or "kendall". Default: "spearman".
#'
#' @return A gene regulatory network represented as an edge list.
#' @export
#' @rdname grn_aracne
#' @importFrom minet aracne
grn_aracne <- function(exp, estimator = "spearman") {
    mi_mat <- minet::build.mim(t(exp), estimator = estimator)
    grn <- minet::aracne(mi_mat, eps = 0.1)
    grn_edges <- cormat_to_edgelist(grn)
    grn_edges <- grn_edges[order(-grn_edges[,3]), ]
    return(grn_edges)
}


#' Infer gene regulatory network with GENIE3
#'
#' @param exp Expression matrix with gene IDs as row names and samples as column names.
#' @param regulators A character vector of regulators (e.g., transcription factors or miRNAs). All regulators must be included in `exp`.
#' @param ... Additional arguments passed to `GENIE3::GENIE3()`.
#'
#' @return A gene regulatory network represented as an edge list.
#' @importFrom GENIE3 GENIE3 getLinkList
#' @rdname grn_genie3
#' @export
#' @importFrom GENIE3 GENIE3 getLinkList
grn_genie3 <- function(exp, regulators = NULL, ...) {
    grn <- GENIE3::GENIE3(as.matrix(exp), regulators = regulators, ...)
    edges <- GENIE3::getLinkList(grn)
    edges[, 1] <- as.character(edges[, 1])
    edges[, 2] <- as.character(edges[, 2])
    edges[, 3] <- as.numeric(edges[, 3])
    colnames(edges) <- c("Node1", "Node2", "Weight")
    edges <- edges[!duplicated(cbind(pmin(edges$Node1, edges$Node2),
                                     pmax(edges$Node1, edges$Node2))), ]
    edges <- edges[order(-edges[,3]), ]
    return(edges)
}


#' Infer gene regulatory network with multiple algorithms and combine results in a list
#'
#' @param exp Expression matrix with gene IDs as row names and samples as column names.
#' @param regulators A character vector of regulators (e.g., transcription factors or miRNAs). All regulators must be included in `exp`.
#' @param ... Additional arguments passed to `GENIE3::GENIE3()`.
#'
#' @return A list of data frames representing edge lists. Each list element is an edge list for a specific method.
#' @rdname grn_combined
#' @export
grn_combined <- function(exp, regulators = NULL, ...) {
    genie3 <- grn_genie3(exp, regulators, ...)
    aracne <- grn_aracne(exp)
    clr <- grn_clr(exp)

    res_list <- list(genie3, aracne, clr)
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






















