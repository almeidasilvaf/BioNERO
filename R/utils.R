
#' Transform a correlation matrix to an edge list
#'
#' @param matrix Symmetrical correlation matrix.
#'
#' @return A 2-column data frame containing node 1, node 2 and edge weight
#' @export
#' @rdname cormat_to_edgelist
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
#' @param edgelist Edge list as a 2-column data frame containing node 1, node 2 and edge weight.
#' @param net_type Type of biological network. One of "gcn", "grn", or "ppi". Default: gcn.
#'
#' @rdname check_sft
#' @export
#' @importFrom igraph graph_from_data_frame as_adjacency_matrix fit_power_law
check_sft <- function(edgelist, net_type = "gcn") {

    # Calculate degree of the resulting graph
    if(net_type == "gcn") {
        graph <- igraph::graph_from_data_frame(edgelist, directed=FALSE)
        adj <- igraph::as_adjacency_matrix(graph, sparse = FALSE)
        diag(adj) <- 0
        degree <- apply(adj, 1, sum, na.rm=TRUE)
    } else if(net_type == "grn") {
        graph <- igraph::graph_from_data_frame(edgelist, directed=TRUE)
        degree <- igraph::degree(graph, mode = "out")
    } else if(net_type == "ppi") {
        graph <- igraph::graph_from_data_frame(edgelist, directed=FALSE)
        degree <- igraph::degree(graph, mode)
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



























