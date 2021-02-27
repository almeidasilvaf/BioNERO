
#' Calculate network statistics
#'
#' @param adj_matrix Adjacency matrix that represents the network.
#' @param net_type One of "gcn" (gene coexpression network), "ppi" (protein-protein interaction), or "grn" (gene regulatory network).
#' @param calculate_additional Logical indicating whether to calculate additional network statistics (betweenness and closeness). Default is FALSE, because it might be time-consuming.
#'
#' @return A list containing the following elements: \itemize{
#'   \item Connectivity
#'   \item ScaledConnectivity
#'   \item ClusterCoef
#'   \item MAR (for gcn only)
#'   \item Density
#'   \item Centralization
#'   \item Heterogeneity (gcn only)
#'   \item nCliques
#'   \item Diameter
#'   \item Betweenness
#'   \item Closeness
#' }
#'
#' @seealso
#'  \code{\link[igraph]{graph_from_adjacency_matrix}},\code{\link[igraph]{cliques}},\code{\link[igraph]{diameter}},\code{\link[igraph]{estimate_betweenness}},\code{\link[igraph]{V}},\code{\link[igraph]{closeness}},\code{\link[igraph]{degree}},\code{\link[igraph]{transitivity}},\code{\link[igraph]{edge_density}},\code{\link[igraph]{centr_degree}}
#'  \code{\link[WGCNA]{fundamentalNetworkConcepts}}
#' @rdname net_stats
#' @export
#' @importFrom igraph graph_from_adjacency_matrix cliques diameter betweenness V closeness degree transitivity edge_density centralization.degree
#' @importFrom WGCNA fundamentalNetworkConcepts
#' @examples
#' \donttest{
#' data(filt.se)
#' set.seed(12)
#' filt.se <- filter_by_variance(filt.se, n=100)
#' gcn <- exp2gcn(filt.se, SFTpower = 19, cor_method = "pearson",
#'                reportPDF = FALSE, net_type = "signed hybrid")
#' stats <- net_stats(gcn$adjacency_matrix, net_type="gcn")
#' }
net_stats <- function(adj_matrix = NULL, net_type=c("gcn", "ppi", "grn"),
                      calculate_additional = FALSE) {

  # For gene coexpression networks
  if(net_type == "gcn") {
    # Create a graph object representing an undirected weighted network
    graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode="undirected",
                                                 diag=FALSE, weighted=TRUE)
    # Calculate stats
    stats <- WGCNA::fundamentalNetworkConcepts(adj_matrix)
    stats$nCliques <- length(igraph::cliques(graph, min=3)) # number of cliques
    stats$diameter <- igraph::diameter(graph, directed=FALSE) # net diameter

    if(calculate_additional) {
      # Calculate vertex betweenness
      betweenness <- igraph::betweenness(graph, directed = FALSE)
      names(betweenness) <- igraph::V(graph)
      stats$betweenness <- as.data.frame(betweenness)

      # Calculate closeness
      closeness <- igraph::closeness(graph, mode="all")
      names(closeness) <- igraph::V(graph)
      stats$closeness <- as.data.frame(closeness)
    }
    # For protein-protein interaction networks
  } else if(net_type == "ppi") {
    # Create a graph object representing an undirected unweighted network
    graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode="undirected",
                                                 diag=FALSE, weighted=FALSE)

    # Calculate stats
    stats <- WGCNA::fundamentalNetworkConcepts(adj_matrix)
    stats$MAR <- NULL # Remove MAR, since it doesn't make sense here
    stats$nCliques <- length(igraph::cliques(graph, min=3)) # number of cliques
    stats$diameter <- igraph::diameter(graph, directed=FALSE) # net diameter

    if(calculate_additional == TRUE) {
      # Calculate vertex betweenness
      betweenness <- igraph::betweenness(graph, directed = FALSE)
      names(betweenness) <- igraph::V(graph)
      stats$betweenness <- as.data.frame(betweenness)

      # Calculate closeness
      closeness <- igraph::closeness(graph, mode="all")
      names(closeness) <- igraph::V(graph)
      stats$closeness <- as.data.frame(closeness)
    }
    # For gene regulatory networks
  } else if(net_type == "grn") {
    graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode="directed",
                                                 diag=FALSE, weighted=FALSE)
    # Calculate stats
    degree <- igraph::degree(graph, mode = "all")
    clustercoef <- igraph::transitivity(graph, mode="global")
    density <- igraph::edge_density(graph)
    centralization <- igraph::centralization.degree(graph)$centralization
    diameter <- igraph::diameter(graph, directed=TRUE)

    if(calculate_additional) {
      # Calculate vertex betweenness
      betweenness <- igraph::betweenness(graph, directed = TRUE)
      names(betweenness) <- igraph::V(graph)
      stats$betweenness <- as.data.frame(betweenness)

      # Calculate closeness
      closeness <- igraph::closeness(graph, mode="all")
      names(closeness) <- igraph::V(graph)
      stats$closeness <- as.data.frame(closeness)

      stats <- list(
        Connectivity = degree,
        ScaledConnectiviy = degree / max(degree),
        ClusterCoef = clustercoef,
        Density = density,
        Centralization = centralization,
        Diameter = diameter,
        Betweenness = betweenness,
        Closeness = closeness)
    } else{
      stats <- list(
        Connectivity = degree,
        ScaledConnectiviy = degree / max(degree),
        ClusterCoef = clustercoef,
        Density = density,
        Centralization = centralization,
        Diameter = diameter)
      }
  } else{
    stop("Please, specify a valid net_type. One of 'ppi', 'gcn' or 'grn'.")
  }
  return(stats)
}











