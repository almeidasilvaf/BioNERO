#' Transform a correlation matrix to an edge list
#'
#' @param matrix Symmetrical correlation matrix.
#'
#' @return A 2-column data frame containing node 1, node 2 and edge weight
#' @export
cormat_to_edgelist <- function(matrix) {
    edgelist <- matrix
    edgelist[lower.tri(edgelist, diag=TRUE)] <- NA
    edgelist <- na.omit(data.frame(as.table(edgelist), stringsAsFactors=FALSE))
    colnames(edgelist) <- c("Node1", "Node2", "Weight")
    return(edgelist)
}


#' Infer gene regulatory network with the Context Likelihood of Relatedness (CLR) algorithm
#'
#' @param exp Expression matrix with gene IDs as row names and samples as column names.
#'
#' @return A gene regulatory network represented as an edge list.
#' @export
#' @importFrom minet build.mim clr
grn_clr <- function(exp) {
    mi_mat <- minet::build.mim(t(exp), estimator = "pearson")
    grn <- minet::clr(mi_mat)
    grn_edges <- cormat_to_edgelist(grn)
    grn_edges <- grn_edges[grn_edges$Weight != 0, ]
    return(grn_edges)
}


#' Infer gene regulatory network with the ARACNE algorithm
#'
#' @param exp Expression matrix with gene IDs as row names and samples as column names.
#'
#' @return A gene regulatory network represented as an edge list.
#' @export
#' @importFrom minet aracne
grn_aracne <- function(exp) {
    mi_mat <- minet::build.mim(t(exp), estimator = "pearson")
    grn <- minet::aracne(mi_mat, eps = 0.1)
    grn_edges <- cormat_to_edgelist(grn)
    grn_edges <- grn_edges[grn_edges$Weight != 0, ]
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
grn_genie3 <- function(exp, regulators = NULL, ...) {
    grn <- GENIE3::GENIE3(as.matrix(exp), regulators = regulators)
    edges <- GENIE3::getLinkList(grn)
    return(edges)
}






















