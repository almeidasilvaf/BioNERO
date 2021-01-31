#' Detect communities in a network
#'
#' @param edgelist Data frame containing the network as an edge list. First column must be node 1 and second column must be node 2. Additional columns will be interpreted as edge attributes and will be modified by this function.
#' @param method Community detection algorithm to be used. Available methods are "infomap", "edge_betweenness", "fast_greedy", "walktrap", "spinglass", "leading_eigen", "louvain", and "label_prop". Default is "infomap".
#' @return A data frame containing node names in the first column, and communities to which nodes belong in the second column.
#'
#' @seealso
#'  \code{\link[igraph]{as_data_frame}},\code{\link[igraph]{simplify}},\code{\link[igraph]{cluster_infomap}},\code{\link[igraph]{cluster_edge_betweenness}},\code{\link[igraph]{cluster_fast_greedy}},\code{\link[igraph]{cluster_walktrap}},\code{\link[igraph]{cluster_spinglass}},\code{\link[igraph]{cluster_leading_eigen}},\code{\link[igraph]{cluster_louvain}},\code{\link[igraph]{cluster_label_prop}}
#' @rdname detect_communities
#' @author Fabricio Almeida-Silva
#' @export
#' @importFrom igraph graph.data.frame simplify cluster_infomap cluster_edge_betweenness cluster_fast_greedy cluster_walktrap cluster_spinglass cluster_leading_eigen cluster_louvain cluster_label_prop
detect_communities <- function(edgelist, method = "infomap") {
    graph <- igraph::graph.data.frame(edgelist, directed = FALSE)
    graph <- igraph::simplify(graph)

    if(method == "infomap") {
        com <- igraph::cluster_infomap(graph)
    } else if(method == "edge_betweenness") {
        com <- igraph::cluster_edge_betweenness(graph)
    } else if(method == "fast_greedy") {
        com <- igraph::cluster_fast_greedy(graph)
    } else if(method == "walktrap") {
        com <- igraph::cluster_walktrap(graph)
    } else if(method == "spinglass") {
        com <- igraph::cluster_spinglass(graph)
    } else if(method == "leading_eigen") {
        com <- igraph::cluster_leading_eigen(graph)
    } else if(method == "louvain") {
        com <- igraph::cluster_louvain(graph)
    } else if(method == "label_prop") {
        com <- igraph::cluster_label_prop(graph)
    } else {
        stop("Please, specify a valid community detection algorithm.")
    }

    df_com <- as.data.frame(list(names = com$names, mem = com$membership))
    return(df_com)
}


#' Create a ggnetwork data frame from an igraph object
#'
#' @param graph Object of class igraph.
#' @param layout Network layout. One of "dh", "drl", "gem", "lgl", "fr", "graphopt", "kk" and "mds". Default is "kk".
#' @param arrow.gap Numeric indicating the distance between nodes and arrows. Default is 0.2
#' @seealso
#'  \code{\link[ggnetwork]{ggnetwork}}
#'  \code{\link[igraph]{layout_with_dh}},\code{\link[igraph]{layout_with_drl}},\code{\link[igraph]{layout_with_gem}},\code{\link[igraph]{layout_with_lgl}},\code{\link[igraph]{layout_with_fr}},\code{\link[igraph]{layout_with_graphopt}},\code{\link[igraph]{layout_with_kk}},\code{\link[igraph]{layout_with_mds}}
#' @rdname igraph2ggnetwork
#' @export
#' @author Fabricio Almeida-Silva
#' @importFrom ggnetwork ggnetwork
#' @importFrom igraph with_dh with_drl with_gem with_lgl with_fr with_graphopt with_kk with_mds
igraph2ggnetwork <- function(graph, layout = "kk", arrow.gap = 0.2) {
    if(layout == "dh") {
        ggnet <- ggnetwork::ggnetwork(graph, layout = igraph::with_dh(),  arrow.gap = arrow.gap)
    } else if(layout == "drl") {
        ggnet <- ggnetwork::ggnetwork(graph, layout = igraph::with_drl(), arrow.gap = arrow.gap)
    } else if(layout == "gem") {
        ggnet <- ggnetwork::ggnetwork(graph, layout = igraph::with_gem(), arrow.gap = arrow.gap)
    } else if(layout == "lgl") {
        ggnet <- ggnetwork::ggnetwork(graph, layout = igraph::with_lgl(), arrow.gap = arrow.gap)
    } else if(layout == "fr") {
        ggnet <- ggnetwork::ggnetwork(graph, layout = igraph::with_fr(), arrow.gap = arrow.gap)
    } else if(layout == "graphopt") {
        ggnet <- ggnetwork::ggnetwork(graph, layout = igraph::with_graphopt(), arrow.gap = arrow.gap)
    } else if(layout == "kk") {
        ggnet <- ggnetwork::ggnetwork(graph, layout = igraph::with_kk(), arrow.gap = arrow.gap)
    } else if(layout == "mds") {
        ggnet <- ggnetwork::ggnetwork(graph, layout = igraph::with_mds(), arrow.gap = arrow.gap)
    } else {
        stop("Please, specify a valid layout.")
    }
    return(ggnet)
}


#' Plot protein-protein interaction network from edge list
#'
#' @param edgelist_int Data frame containing the edge list for the PPI network. First column is the protein 1 and second column is the protein 2. All other columns are interpreted as edge attributes.
#' @param detect_communities Logical indicating whether to detect communities or not. Default is TRUE.
#' @param clustering_method Community detection algorithm to be used. Available methods are "infomap", "edge_betweenness", "fast_greedy", "walktrap", "spinglass", "leading_eigen", "louvain", and "label_prop". Default is "infomap".
#' @param show_labels Character indicating which nodes will be labeled. One of "all", "allhubs", "tophubs", or "none".
#' @param top_n_hubs Number of top hubs to be labeled. It is only valid if \code{show_labels} equals "tophubs". Default is 5.
#' @param interactive Logical indicating whether the network should be interactive or not. Default is FALSE.
#' @seealso
#'  \code{\link[igraph]{as_data_frame}},\code{\link[igraph]{degree}},\code{\link[igraph]{simplify}},\code{\link[igraph]{gorder}}
#'  \code{\link[networkD3]{igraph_to_networkD3}},\code{\link[networkD3]{forceNetwork}}
#'  \code{\link[ggnetwork]{geom_edges}},\code{\link[ggnetwork]{geom_nodes}},\code{\link[ggnetwork]{geom_nodetext}},\code{\link[ggnetwork]{theme_blank}},\code{\link[ggnetwork]{geom_nodetext_repel}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes_}},
#' @rdname plot_ppi
#' @author Fabricio Almeida-Silva
#' @export
#' @importFrom igraph graph_from_data_frame degree simplify vcount
#' @importFrom networkD3 igraph_to_networkD3 forceNetwork
#' @importFrom ggnetwork ggnetwork geom_edges geom_nodes geom_nodetext theme_blank geom_nodelabel_repel unit
#' @importFrom ggplot2 ggplot aes_ guides
plot_ppi <- function(edgelist_int, detect_communities = TRUE,
                     clustering_method = "infomap", show_labels = "tophubs",
                     top_n_hubs = 5, interactive = FALSE) {
    requireNamespace("intergraph", quietly=TRUE)
    nod_at <- data.frame(Gene = unique(c(as.character(edgelist_int[,1]), as.character(edgelist_int[,2]))),
                         stringsAsFactors = FALSE)

    # Add communities
    if(detect_communities) {
        clusters <- detect_communities(edgelist_int, method = clustering_method)
        nod_at <- merge(nod_at, clusters, by.x="Gene", by.y="names")
    }

    # Add degree
    g <- igraph::graph_from_data_frame(d = edgelist_int, directed = FALSE)
    g_degree <- as.data.frame(igraph::degree(g), stringsAsFactors=FALSE); colnames(g_degree) <- "Degree"
    nod_at <- merge(nod_at, g_degree, by.x="Gene", by.y="row.names")
    nod_at <- nod_at[order(-nod_at$Degree), ]

    # Add hub gene status
    hubs <- nod_at[order(nod_at$Degree, decreasing = TRUE), ]
    hubs <- hubs[1:(nrow(hubs) / 10), ]
    nod_at$isHub <- ifelse(nod_at$Gene %in% hubs[,1], TRUE, FALSE)

    # Should the network be interactive?
    if(interactive) {
        graph <- igraph::simplify(igraph::graph_from_data_frame(d = edgelist_int, vertices = nod_at, directed=FALSE))

        graph_d3 <- networkD3::igraph_to_networkD3(graph, group = nod_at$mem)
        graph_d3$nodes <- merge(graph_d3$nodes, nod_at, by.x="name", by.y="Gene", sort = FALSE)
        p <- networkD3::forceNetwork(Links = graph_d3$links, Nodes = graph_d3$nodes,
                                     Source = 'source', Target = 'target',
                                     NodeID = 'name', Group = 'group',
                                     Nodesize = 'Degree', height=900, width=1200,
                                     opacity=0.8, zoom = TRUE, fontSize = 20)

    } else { #Static network
        # Define plotting parameters
        if(show_labels == "all") {
            graph <- igraph::graph_from_data_frame(d = edgelist_int, vertices = nod_at, directed = FALSE)
            nvertices <- igraph::vcount(graph)
            if(nvertices > 200) {
                message("WARNING: Graph has more than 400 vertices. Consider displaying labels of hubs or top hubs only.
                  Number of vertices:", nvertices)
            }
            graph <- igraph::simplify(graph)
            n <- ggnetwork::ggnetwork(graph)
            n$Cluster <- as.factor(n$mem)
            p <- ggplot2::ggplot(n, ggplot2::aes_(x = ~x, y = ~y, xend = ~xend, yend = ~yend)) +
                ggnetwork::geom_edges(color = "grey75", alpha = 0.5, show.legend = FALSE) +
                ggnetwork::geom_nodes(ggplot2::aes_(size = ~Degree, color = ~Cluster)) +
                ggplot2::guides(color = FALSE) +
                ggnetwork::geom_nodetext(ggplot2::aes_(label = ~name, size = 0.4 * ~Degree)) +
                ggnetwork::theme_blank()

        } else if(show_labels == "allhubs") {
            graph <- igraph::graph_from_data_frame(d = edgelist_int, vertices = nod_at, directed = FALSE)
            graph <- igraph::simplify(graph)
            n <- ggnetwork::ggnetwork(graph)
            n$Cluster <- as.factor(n$mem)
            p <- ggplot2::ggplot(n, ggplot2::aes_(x = ~x, y = ~y, xend = ~xend, yend = ~yend)) +
                ggnetwork::geom_edges(color = "grey75", alpha = 0.5, show.legend = FALSE) +
                ggnetwork::geom_nodes(ggplot2::aes_(size = ~Degree, color = ~Cluster)) +
                ggplot2::guides(color = FALSE) +
                ggnetwork::geom_nodelabel_repel(ggplot2::aes_(label = ~name, color = ~isHub),
                                                box.padding = ggnetwork::unit(1, "lines"),
                                                data = function(x) { x[ x$isHub, ]}, show.legend = FALSE) +
                ggnetwork::theme_blank()
        } else if(show_labels == "tophubs") {
            tophubs <- nod_at[nod_at$isHub == TRUE, 1][1:top_n_hubs]
            nod_at$isTopHub <- ifelse(nod_at$Gene %in% tophubs, TRUE, FALSE)
            graph <- igraph::graph_from_data_frame(d = edgelist_int, vertices = nod_at, directed = FALSE)
            graph <- igraph::simplify(graph)
            n <- ggnetwork::ggnetwork(graph)
            n$Cluster <- as.factor(n$mem)
            p <- ggplot2::ggplot(n, ggplot2::aes_(x = ~x, y = ~y, xend = ~xend, yend = ~yend)) +
                ggnetwork::geom_edges(color = "grey75", alpha = 0.5, show.legend = FALSE) +
                ggnetwork::geom_nodes(ggplot2::aes_(size = ~Degree, color = ~Cluster)) +
                ggplot2::guides(color = FALSE) +
                ggnetwork::geom_nodelabel_repel(ggplot2::aes_(label = ~name, color = ~isHub),
                                                box.padding = ggnetwork::unit(1, "lines"),
                                                data = function(x) { x[ x$isTopHub, ]}, show.legend = FALSE) +
                ggnetwork::theme_blank()
        } else if(show_labels == "none") {
            graph <- igraph::graph_from_data_frame(d = edgelist_int, vertices = nod_at, directed = FALSE)
            graph <- igraph::simplify(graph)
            n <- ggnetwork::ggnetwork(graph)
            n$Cluster <- as.factor(n$mem)
            p <- ggplot2::ggplot(n, ggplot2::aes_(x = ~x, y = ~y, xend = ~xend, yend = ~yend)) +
                ggnetwork::geom_edges(color = "grey75", alpha = 0.5, show.legend = FALSE) +
                ggnetwork::geom_nodes(ggplot2::aes_(size = ~Degree, color = ~Cluster)) +
                ggplot2::guides(color = FALSE) +
                ggnetwork::theme_blank()
        }
    }

    return(p)
}



#' Plot gene regulatory network from edge list
#'
#' @param edgelist_grn Data frame containing the edge list for the GRN network. First column is the TF and second column is the target gene. All other columns are interpreted as edge attributes.
#' @param show_labels Character indicating which nodes will be labeled. One of "all", "allhubs", "tophubs", or "none".
#' @param top_n_hubs Number of top hubs to be labeled. It is only valid if \code{show_labels} equals "tophubs". Default is 5.
#' @param interactive Logical indicating whether the network should be interactive or not. Default is FALSE.
#' @param layout Network layout. One of "dh", "drl", "gem", "lgl", "fr", "graphopt", "kk" and "mds". Default is "kk".
#' @param arrow.gap Numeric indicating the distance between nodes and arrows. Default is 0.2.
#'
#' @return A ggplot object containing the network.
#' @seealso
#'  \code{\link[igraph]{as_data_frame}},\code{\link[igraph]{degree}},\code{\link[igraph]{gorder}}
#'  \code{\link[networkD3]{igraph_to_networkD3}},\code{\link[networkD3]{forceNetwork}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes_}},\code{\link[ggplot2]{reexports}},\code{\link[ggplot2]{scale_manual}}
#'  \code{\link[ggnetwork]{geom_edges}}, \code{\link[ggnetwork]{geom_nodes}},\code{\link[ggnetwork]{geom_nodetext}},\code{\link[ggnetwork]{theme_blank}},\code{\link[ggnetwork]{geom_nodetext_repel}}
#'  \code{\link[ggnewscale]{new_scale}}
#' @rdname plot_grn
#' @author Fabricio Almeida-Silva
#' @export
#' @importFrom igraph graph_from_data_frame degree vcount
#' @importFrom networkD3 igraph_to_networkD3 forceNetwork
#' @importFrom ggplot2 ggplot aes_ arrow scale_color_manual guides
#' @importFrom ggnetwork geom_edges unit geom_nodes geom_nodetext theme_blank geom_nodelabel_repel
#' @importFrom ggnewscale new_scale_color
plot_grn <- function(edgelist_grn, show_labels = "tophubs", top_n_hubs = 5,
                     interactive = FALSE, layout = "kk", arrow.gap = 0.01) {
    requireNamespace("intergraph", quietly=TRUE)

    if(ncol(edgelist_grn) == 3) {
        colnames(edgelist_grn)[3] <- "Regulation"
        palette <- c("red", "blue", "grey75")
        palette2 <- c("gold2", "forestgreen", "grey75")
        showlegend <- TRUE
    } else {
        edgelist_grn$Regulation <- "none"
        palette <- c("grey75", "grey76", "grey77")
        palette2 <- c("gold2", "forestgreen")
        showlegend <- FALSE
    }

    # Start data frame of node attributes
    nod_at <- data.frame(Gene = unique(c(as.character(edgelist_grn[,1]), as.character(edgelist_grn[,2]))),
                         stringsAsFactors = FALSE)

    # Add degree
    g <- igraph::graph_from_data_frame(d = edgelist_grn, directed = TRUE)
    g_degree <- as.data.frame(igraph::degree(g, mode = "out"), stringsAsFactors=FALSE); colnames(g_degree) <- "Degree"
    nod_at <- merge(nod_at, g_degree, by.x="Gene", by.y="row.names")
    nod_at <- nod_at[order(-nod_at$Degree), ]

    # Add hub gene status
    hubs <- nod_at[order(nod_at$Degree, decreasing = TRUE), ]
    hubs <- hubs[1:(nrow(hubs) / 10), ]
    nod_at$isHub <- ifelse(nod_at$Gene %in% hubs[,1], TRUE, FALSE)

    # Add classification for TF and target gene
    nod_at$Molecule <- ifelse(nod_at$Gene %in% edgelist_grn[,1], "TF", "target")

    # Should the network be interactive?
    if(interactive) {
        graph <- igraph::graph_from_data_frame(d = edgelist_grn, vertices = nod_at, directed=TRUE)

        graph_d3 <- networkD3::igraph_to_networkD3(graph, group = nod_at$Molecule)
        graph_d3$nodes <- merge(graph_d3$nodes, nod_at, by.x="name", by.y="Gene", sort = FALSE)
        my_color <- 'd3.scaleOrdinal() .domain(["TF", "target"]) .range(["forestgreen", "orange"])'
        p <- networkD3::forceNetwork(Links = graph_d3$links, Nodes = graph_d3$nodes,
                                     Source = 'source', Target = 'target',
                                     NodeID = 'name', Group = 'group',
                                     Value = 'value',
                                     linkColour = ifelse(graph_d3$links$value == "positive", "red", "blue"),
                                     colourScale = my_color,
                                     Nodesize = 'Degree', height=900, width=1200,
                                     opacity=1, zoom = TRUE, fontSize = 20, legend=TRUE)

    } else { #Static network
        # Define plotting parameters
        if(show_labels == "all") {
            graph <- igraph::graph_from_data_frame(d = edgelist_grn, vertices = nod_at, directed = TRUE)
            nvertices <- igraph::vcount(graph)
            if(nvertices > 200) {
                message("WARNING: Graph has more than 400 vertices. Consider displaying labels of hubs or top hubs only.
                  Number of vertices:", nvertices)
            }
            n <- igraph2ggnetwork(graph, layout = layout, arrow.gap = arrow.gap)

            p <- ggplot2::ggplot(n, ggplot2::aes_(x = ~x, y = ~y, xend = ~xend, yend = ~yend)) +
                ggnetwork::geom_edges(ggplot2::aes_(color = ~Regulation), alpha = 0.5,
                                      arrow = ggplot2::arrow(length = ggnetwork::unit(0.1, "lines"), type = "closed"),
                                      curvature = 0.1, show.legend = showlegend) +
                ggplot2::scale_color_manual(values = palette) +
                ggnewscale::new_scale_color() +
                ggnetwork::geom_nodes(ggplot2::aes_(color = ~Molecule, size = ~Degree, shape = ~Molecule)) +
                ggplot2::scale_color_manual(values = palette2) +
                ggnetwork::geom_nodetext(ggplot2::aes_(label = ~name, size = 0.5 * ~Degree), vjust = -1) +
                ggnetwork::theme_blank()

        } else if(show_labels == "allhubs") {
            graph <- igraph::graph_from_data_frame(d = edgelist_grn, vertices = nod_at, directed = TRUE)
            n <- igraph2ggnetwork(graph, layout = layout, arrow.gap = arrow.gap)
            p <- ggplot2::ggplot(n, ggplot2::aes_(x = ~x, y = ~y, xend = ~xend, yend = ~yend)) +
                ggnetwork::geom_edges(ggplot2::aes_(color = ~Regulation), alpha = 0.5,
                                      arrow = ggplot2::arrow(length = ggnetwork::unit(0.1, "lines"), type = "closed"),
                                      curvature = 0.1, show.legend = showlegend) +
                ggplot2::scale_color_manual(values = palette) +
                ggnewscale::new_scale_color() +
                ggnetwork::geom_nodes(ggplot2::aes_(color = ~Molecule, size = ~Degree, shape = ~Molecule)) +
                ggplot2::scale_color_manual(values = palette2) +
                ggnewscale::new_scale_color() +
                ggnetwork::geom_nodelabel_repel(ggplot2::aes_(label = ~name, color = ~isHub),
                                                box.padding = ggnetwork::unit(1, "lines"),
                                                data = function(x) { x[ x$isHub, ]}, show.legend = FALSE) +
                ggnetwork::theme_blank()
        } else if(show_labels == "tophubs") {
            tophubs <- nod_at[nod_at$isHub == TRUE, 1][1:top_n_hubs]
            nod_at$isTopHub <- ifelse(nod_at$Gene %in% tophubs, TRUE, FALSE)
            graph <- igraph::graph_from_data_frame(d = edgelist_grn, vertices = nod_at, directed = TRUE)
            n <- igraph2ggnetwork(graph, layout = layout, arrow.gap = arrow.gap)
            p <- ggplot2::ggplot(n, ggplot2::aes_(x = ~x, y = ~y, xend = ~xend, yend = ~yend)) +
                ggnetwork::geom_edges(ggplot2::aes_(color = ~Regulation), alpha = 0.5,
                                      arrow = ggplot2::arrow(length = ggnetwork::unit(0.1, "lines"), type = "closed"),
                                      curvature = 0.1, show.legend = showlegend) +
                ggplot2::scale_color_manual(values = palette) +
                ggnewscale::new_scale_color() +
                ggnetwork::geom_nodes(ggplot2::aes_(color = ~Molecule, size = ~Degree, shape = ~Molecule)) +
                ggplot2::scale_color_manual(values = palette2) +
                ggnewscale::new_scale_color() +
                ggplot2::guides(color = FALSE) +
                ggnetwork::geom_nodelabel_repel(ggplot2::aes_(label = ~name, color = ~isHub),
                                                box.padding = ggnetwork::unit(1, "lines"),
                                                data = function(x) { x[ x$isTopHub, ]}, show.legend = FALSE) +
                ggnetwork::theme_blank()
        } else if(show_labels == "none") {
            graph <- igraph::graph_from_data_frame(d = edgelist_grn, vertices = nod_at, directed = TRUE)
            n <- igraph2ggnetwork(graph, layout = layout, arrow.gap = arrow.gap)
            p <- ggplot2::ggplot(n, ggplot2::aes_(x = ~x, y = ~y, xend = ~xend, yend = ~yend)) +
                ggnetwork::geom_edges(ggplot2::aes_(color = ~Regulation), alpha = 0.5,
                                      arrow = ggplot2::arrow(length = ggnetwork::unit(0.1, "lines"), type = "closed"),
                                      curvature = 0.1, show.legend = showlegend) +
                ggplot2::scale_color_manual(values = palette) +
                ggnewscale::new_scale_color() +
                ggnetwork::geom_nodes(ggplot2::aes_(color = ~Molecule, size = ~Degree, shape = ~Molecule)) +
                ggplot2::scale_color_manual(values = palette2) +
                ggnetwork::theme_blank()
        }
    }

    return(p)
}


#' Plot gene coexpression network from edge list
#'
#' @param edgelist_gcn Data frame containing the edge list for the GCN. The edge list can be generated with \code{get_edge_list()}.
#' @param net List object returned by \code{exp2net}.
#' @param color_by How should nodes be colored? It must be either "module" (nodes will have the colors of their modules) or a 2-column data frame containing genes in the first column and a custom gene annotation in the second column. Default: "module".
#' @param hubs Data frame containing hub genes in the first column, their modules in the second column, and intramodular connectivity in the third column.
#' @param show_labels Character indicating which nodes will be labeled. One of "all", "allhubs", "tophubs", or "none". Default: tophubs.
#' @param top_n_hubs Number of top hubs to be labeled. It is only valid if \code{show_labels} equals "tophubs". Default is 5.
#' @param interactive Logical indicating whether the network should be interactive or not. Default is FALSE.
#' @seealso
#'  \code{\link[igraph]{simplify}},\code{\link[igraph]{as_data_frame}},\code{\link[igraph]{gorder}}
#'  \code{\link[networkD3]{igraph_to_networkD3}},\code{\link[networkD3]{forceNetwork}}
#'  \code{\link[ggnetwork]{geom_edges}},\code{\link[ggnetwork]{geom_nodes}},\code{\link[ggnetwork]{geom_nodetext}},\code{\link[ggnetwork]{theme_blank}},\code{\link[ggnetwork]{geom_nodetext_repel}}
#'  \code{\link[ggplot2]{ggplot}},\code{\link[ggplot2]{aes_}}
#' @rdname plot_gcn
#' @author Fabricio Almeida-Silva
#' @export
#' @importFrom igraph simplify graph_from_data_frame
#' @importFrom networkD3 igraph_to_networkD3 forceNetwork
#' @importFrom ggnetwork ggnetwork geom_edges geom_nodes geom_nodetext theme_blank geom_nodelabel_repel unit
#' @importFrom ggplot2 ggplot aes_ guides
plot_gcn <- function(edgelist_gcn, net, color_by="module", hubs = NULL,
                     show_labels = "tophubs", top_n_hubs = 5,
                     interactive = FALSE) {
    requireNamespace("intergraph", quietly=TRUE)

    if(is.null(hubs) | is.null(edgelist_gcn)) {
        stop("Arguments 'edgelist_gcn' and 'hubs' are mandatory for this network.")
    }

    # How should genes be colored?
    if(is.data.frame(color_by)) {
        gene_annotation <- color_by # color by custom annotation
    } else {
        gene_annotation <- net$genes_and_modules # color by module color
    }
    kIN <- net$kIN

    # Create a data frame of nodes and node attributes
    geneIDs <- unique(c(as.character(edgelist_gcn[,1]), as.character(edgelist_gcn[,2])))
    nod_at <- data.frame(Gene = geneIDs, stringsAsFactors = FALSE)
    nod_at$Class <- as.factor(gene_annotation[gene_annotation[,1] %in% nod_at$Gene, 2])
    nod_at$Degree <- kIN$kWithin[rownames(kIN) %in% nod_at$Gene]
    nod_at$isHub <- ifelse(nod_at$Gene %in% hubs[,1], TRUE, FALSE)
    nod_at <- nod_at[order(nod_at$Class, -nod_at$Degree), ]

    # Should the network be interactive?
    if(interactive) {
        graph <- igraph::simplify(igraph::graph_from_data_frame(d = edgelist_gcn, vertices = nod_at, directed=FALSE))
        graph_d3 <- networkD3::igraph_to_networkD3(graph, group = nod_at$Class)
        graph_d3$nodes <- merge(graph_d3$nodes, nod_at, by.x="name", by.y="Gene", sort = FALSE)
        p <- networkD3::forceNetwork(Links = graph_d3$links, Nodes = graph_d3$nodes,
                                     Source = 'source', Target = 'target',
                                     NodeID = 'name', Group = 'group',
                                     Nodesize = 'Degree', height=900, width=1200,
                                     opacity=0.8, zoom = TRUE, fontSize = 13)
    } else {
        # Define plotting parameters
        # Handle gene coloring based on number of annotation classes
        if(nlevels(nod_at$Class) == 1) {
            nodecol <- levels(nod_at$Class)
            add_nodes <- ggnetwork::geom_nodes(ggplot2::aes_(size = ~Degree, alpha = ~Degree),
                                               color = nodecol)
        } else {
            add_nodes <- ggnetwork::geom_nodes(ggplot2::aes_(size = ~Degree, alpha = ~Degree,
                                                             color = ~Class))
        }

        if(is.data.frame(color_by)) {
            palette <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF",
                         "#9467BDFF", "#8C564BFF", "#E377C2FF", "#7F7F7FFF",
                         "#BCBD22FF", "#17BECFFF", "#AEC7E8FF", "#FFBB78FF",
                         "#98DF8AFF", "#FF9896FF", "#C5B0D5FF", "#C49C94FF",
                         "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF")
            scale_color <- ggplot2::scale_color_manual(values = palette[1:nlevels(nod_at$Class)])
        } else {
            scale_color <- ggplot2::scale_color_manual(values = levels(nod_at$Class))
        }

        # Handle node labeling based on which labels to display (top hubs, hubs or all genes)
        if(show_labels == "all") {
            nod_at$Degree2 <- nod_at$Degree * 0.4
            add_nodelabel <- ggnetwork::geom_nodetext(ggplot2::aes_(label = ~name, size = ~Degree2),
                                                      show.legend = FALSE)
        } else if(show_labels == "allhubs") {
            add_nodelabel <- ggnetwork::geom_nodelabel_repel(ggplot2::aes_(label = ~name),
                                                             color="azure4",
                                                             box.padding = ggnetwork::unit(1, "lines"),
                                                             data = function(x) { x[ x$isHub, ]},
                                                             show.legend = FALSE, max.overlaps=Inf)
        } else if(show_labels == "tophubs") {
            tophubs <- nod_at[nod_at$isHub == TRUE, 1][1:top_n_hubs]
            nod_at$isTopHub <- ifelse(nod_at$Gene %in% tophubs, TRUE, FALSE)
            add_nodelabel <- ggnetwork::geom_nodelabel_repel(ggplot2::aes_(label = ~name),
                                                             color="azure4",
                                                             box.padding = ggnetwork::unit(1, "lines"),
                                                             data = function(x) { x[ x$isTopHub, ]},
                                                             show.legend = FALSE)

        } else if(show_labels == "none") {
            add_nodelabel <- NULL
        } else {
            stop("Please, specify a valid option for 'show_labels'.")
        }

        # Create graph object
        graph <- igraph::graph_from_data_frame(d=edgelist_gcn, vertices=nod_at, directed=FALSE)
        n <- igraph2ggnetwork(graph, arrow.gap=0)

        # Plot graph
        p <- ggplot2::ggplot(n, ggplot2::aes_(x = ~x, y = ~y, xend = ~xend, yend = ~yend)) +
            ggnetwork::geom_edges(color = "grey75", alpha = 0.5, show.legend=FALSE) +
            add_nodes +
            scale_color +
            add_nodelabel +
            ggnetwork::theme_blank()
    }

    return(p)
}
