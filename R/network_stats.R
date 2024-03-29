
#' Calculate network statistics
#'
#' @param adj_matrix Adjacency matrix that represents the network.
#' @param net_type One of "gcn" (gene coexpression network),
#' "ppi" (protein-protein interaction), or "grn" (gene regulatory network).
#' @param calculate_additional Logical indicating whether to calculate
#' additional network statistics (betweenness and closeness).
#' Default is FALSE.
#'
#' @return A list containing the following elements: \itemize{
#'   \item Connectivity
#'   \item ScaledConnectivity
#'   \item ClusterCoef
#'   \item MAR (for gcn only)
#'   \item Density
#'   \item Centralization
#'   \item Heterogeneity (gcn only)
#'   \item Diameter
#'   \item Betweenness
#'   \item Closeness
#' }
#'
#' @seealso
#'  \code{\link[igraph]{graph_from_adjacency_matrix}},
#'  \code{\link[igraph]{cliques}},\code{\link[igraph]{diameter}},
#'  \code{\link[igraph]{estimate_betweenness}},\code{\link[igraph]{V}},
#'  \code{\link[igraph]{closeness}},\code{\link[igraph]{degree}},
#'  \code{\link[igraph]{transitivity}},\code{\link[igraph]{edge_density}},
#'  \code{\link[igraph]{centr_degree}}
#'  \code{\link[WGCNA]{fundamentalNetworkConcepts}}
#' @rdname net_stats
#' @export
#' @importFrom igraph graph_from_adjacency_matrix cliques diameter
#' betweenness V closeness degree transitivity edge_density
#' centr_degree
#' @importFrom WGCNA fundamentalNetworkConcepts
#' @examples
#' \donttest{
#' data(filt.se)
#' set.seed(12)
#' filt.se <- exp_preprocess(
#'     filt.se, Zk_filtering = FALSE, variance_filter = TRUE, n = 200
#' )
#' gcn <- exp2gcn(
#'     filt.se, SFTpower = 7, cor_method = "pearson", net_type = "signed hybrid"
#' )
#' stats <- net_stats(gcn$adjacency_matrix, net_type = "gcn")
#' }
net_stats <- function(adj_matrix = NULL, net_type = c("gcn", "ppi", "grn"),
                      calculate_additional = FALSE) {

    weighted <- NULL
    mode <- "undirected"
    directed <- FALSE
    if(net_type == "gcn") {
        weighted <- TRUE
    } else if(net_type == "grn") {
        mode <- "directed"
        directed <- TRUE
    } else if(net_type == "ppi") {
        a <- NULL
    } else {
        stop("Please, specify a valid network type. One of 'gcn', 'ppi' or 'grn'.")
    }
    graph <- igraph::graph_from_adjacency_matrix(
        adj_matrix, mode = mode, diag = FALSE, weighted = weighted
    )
    diam <- igraph::diameter(graph, directed = directed)
    if(net_type != "grn") {
        stats <- WGCNA::fundamentalNetworkConcepts(adj_matrix)
        stats$Diameter <- diam
        if(net_type == "ppi") { stats$MAR <- NULL }
    } else {
        degree <- igraph::degree(graph, mode = "all")
        clustercoef <- igraph::transitivity(graph, type = "global")
        density <- igraph::edge_density(graph)
        centralization <- igraph::centr_degree(graph)$centralization
        stats <- list(
            Connectivity = degree,
            ScaledConnectiviy = degree / max(degree),
            ClusterCoef = clustercoef,
            Density = density,
            Centralization = centralization,
            Diameter = diam
        )
    }
    if(calculate_additional) {
        betweenness <- igraph::betweenness(graph, directed = directed)
        names(betweenness) <- igraph::V(graph)
        closeness <- igraph::closeness(graph, mode="all")
        names(closeness) <- igraph::V(graph)
        stats$betweenness <- as.data.frame(betweenness)
        stats$closeness <- as.data.frame(closeness)
    }
    return(stats)
}











