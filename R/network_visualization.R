#' Detect communities in a network
#'
#' @param edgelist Data frame containing the network as an edge list. First column must be node 1 and second column must be node 2. Additional columns will be interpreted as edge attributes and will be modified by this function.
#' @param method Community detection algorithm to be used. Available methods are "infomap", "edge_betweenness", "fast_greedy", "walktrap", "spinglass", "leading_eigen", "louvain", and "label_prop". Default is "infomap".
#' @param directed Logical indicating whether the network is directed (GRN only) or not (GCN and PPI networks). Default: TRUE.
#' @return A data frame containing node names in the first column, and communities to which nodes belong in the second column.
#'
#' @seealso
#'  \code{\link[igraph]{as_data_frame}},\code{\link[igraph]{simplify}},\code{\link[igraph]{cluster_infomap}},\code{\link[igraph]{cluster_edge_betweenness}},\code{\link[igraph]{cluster_fast_greedy}},\code{\link[igraph]{cluster_walktrap}},\code{\link[igraph]{cluster_spinglass}},\code{\link[igraph]{cluster_leading_eigen}},\code{\link[igraph]{cluster_louvain}},\code{\link[igraph]{cluster_label_prop}}
#' @rdname detect_communities
#' @author Fabricio Almeida-Silva
#' @export
#' @importFrom igraph graph.data.frame simplify cluster_infomap cluster_edge_betweenness cluster_fast_greedy cluster_walktrap cluster_spinglass cluster_leading_eigen cluster_louvain cluster_label_prop
#' @examples
#' data(filt.se)
#' tfs <- sample(rownames(filt.se), size=50, replace=FALSE)
#' grn_edges <- grn_clr(filt.se, regulators = tfs)
#' com <- detect_communities(grn_edges, directed=TRUE)
detect_communities <- function(edgelist, method = "infomap", directed=TRUE) {
    graph <- igraph::graph.data.frame(edgelist, directed = directed)
    graph <- igraph::simplify(graph)

    if(method == "infomap") {
        com <- igraph::cluster_infomap(graph, modularity = FALSE)
    } else if(method == "edge_betweenness") {
        com <- igraph::cluster_edge_betweenness(graph, modularity = FALSE)
    } else if(method == "fast_greedy") {
        com <- igraph::cluster_fast_greedy(graph, modularity = FALSE)
    } else if(method == "walktrap") {
        com <- igraph::cluster_walktrap(graph, modularity = FALSE)
    } else if(method == "spinglass") {
        com <- igraph::cluster_spinglass(graph, modularity = FALSE)
    } else if(method == "leading_eigen") {
        com <- igraph::cluster_leading_eigen(graph, modularity = FALSE)
    } else if(method == "louvain") {
        com <- igraph::cluster_louvain(graph, modularity = FALSE)
    } else if(method == "label_prop") {
        com <- igraph::cluster_label_prop(graph, modularity = FALSE)
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
#' @param arrow.gap Numeric indicating the distance between nodes and arrows. Default is 0.2.
#' @seealso
#'  \code{\link[ggnetwork]{ggnetwork}}
#'  \code{\link[igraph]{layout_with_dh}},\code{\link[igraph]{layout_with_drl}},\code{\link[igraph]{layout_with_gem}},\code{\link[igraph]{layout_with_lgl}},\code{\link[igraph]{layout_with_fr}},\code{\link[igraph]{layout_with_graphopt}},\code{\link[igraph]{layout_with_kk}},\code{\link[igraph]{layout_with_mds}}
#' @rdname igraph2ggnetwork
#' @noRd
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
#' @param color_by How should nodes be colored? It must be either "community" or a 2-column data frame containing proteins in the first column and a custom annotation in the second column. If "community", a clustering algorithm will be applied. Default: "community".
#' @param clustering_method Community detection algorithm to be used. Available methods are "infomap", "edge_betweenness", "fast_greedy", "walktrap", "spinglass", "leading_eigen", "louvain", and "label_prop". Default is "infomap".
#' @param show_labels Character indicating which nodes will be labeled. One of "all", "allhubs", "tophubs", or "none".
#' @param top_n_hubs Number of top hubs to be labeled. It is only valid if \code{show_labels} equals "tophubs". Default is 5.
#' @param interactive Logical indicating whether the network should be interactive or not. Default is FALSE.
#' @param add_color_legend Logical indicating whether to add a color legend for nodes. Default: TRUE.
#' @seealso
#'  \code{\link[igraph]{as_data_frame}},\code{\link[igraph]{degree}},\code{\link[igraph]{simplify}},\code{\link[igraph]{gorder}}
#'  \code{\link[networkD3]{igraph_to_networkD3}},\code{\link[networkD3]{forceNetwork}}
#' @rdname plot_ppi
#' @author Fabricio Almeida-Silva
#' @export
#' @importFrom igraph graph_from_data_frame degree simplify vcount
#' @importFrom networkD3 igraph_to_networkD3 forceNetwork
#' @importFrom ggnetwork ggnetwork geom_edges geom_nodes geom_nodetext theme_blank geom_nodelabel_repel unit
#' @importFrom ggplot2 ggplot aes_ guides
#' @examples
#' ppi_edges <- igraph::get.edgelist(igraph::barabasi.game(n=50, directed=FALSE))
#' p <- plot_ppi(ppi_edges, add_color_legend = FALSE)
plot_ppi <- function(edgelist_int, color_by = "community",
                     clustering_method = "infomap",
                     show_labels = "tophubs",
                     top_n_hubs = 5, interactive = FALSE,
                     add_color_legend=TRUE) {
    requireNamespace("intergraph", quietly=TRUE)
    # How should genes be colored?
    if(is.data.frame(color_by)) {
        prot_annotation <- color_by # color by custom annotation
    } else {
        prot_annotation <- detect_communities(edgelist_int,
                                              method = clustering_method,
                                              directed = FALSE)
    }

    # Get degree of nodes and find hubs
    h <- get_hubs_ppi(edgelist_int, return_degree = TRUE)

    # Create a data frame of nodes and node attributes
    protIDs <- unique(c(as.character(edgelist_int[,1]), as.character(edgelist_int[,2])))
    nod_at <- data.frame(Protein = protIDs, stringsAsFactors = FALSE)
    nod_at$Class <- as.factor(prot_annotation[prot_annotation[,1] %in% nod_at$Protein, 2])
    nod_at$Degree <- h$Degree$Degree[h$Degree$Protein %in% nod_at$Protein]
    nod_at$isHub <- ifelse(nod_at$Protein %in% h$Hubs$Protein, TRUE, FALSE)
    nod_at <- nod_at[order(nod_at$Class, -nod_at$Degree), ]

    # Should the network be interactive?
    if(interactive) {
        graph <- igraph::graph_from_data_frame(d = edgelist_int,
                                               vertices = nod_at, directed=FALSE)
        graph <- igraph::simplify(graph)
        graph_d3 <- networkD3::igraph_to_networkD3(graph, group = nod_at$mem)
        graph_d3$nodes <- merge(graph_d3$nodes, nod_at, by.x="name", by.y="Gene", sort = FALSE)
        p <- networkD3::forceNetwork(Links = graph_d3$links, Nodes = graph_d3$nodes,
                                     Source = 'source', Target = 'target',
                                     NodeID = 'name', Group = 'group',
                                     Nodesize = 'Degree', height=900, width=1200,
                                     opacity=0.8, zoom = TRUE, fontSize = 20)

    } else { #Static network
        # Define plotting parameters
        if(nlevels(nod_at$Class) <= 20) {
            palette <- custom_palette(1)
        } else {
            n <- nlevels(nod_at$Class)
            palette <- colorRampPalette(custom_palette(1))(n)
        }

        # Handle legend
        if(add_color_legend) {
            set_legend <- NULL
        } else {
            set_legend <- ggplot2::guides(color=FALSE)
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
            nod_at$isTopHub <- ifelse(nod_at$Protein %in% tophubs, TRUE, FALSE)
            add_nodelabel <- ggnetwork::geom_nodelabel_repel(ggplot2::aes_(label = ~name),
                                                             color="azure4",
                                                             box.padding = ggnetwork::unit(1, "lines"),
                                                             data = function(x) { x[ x$isTopHub, ]},
                                                             show.legend = FALSE, max.overlaps=Inf)

        } else if(show_labels == "none") {
            add_nodelabel <- NULL
        } else {
            stop("Please, specify a valid option for 'show_labels'.")
        }

        # Create graph object
        graph <- igraph::simplify(igraph::graph_from_data_frame(d=edgelist_int,
                                               vertices=nod_at, directed=FALSE))
        n <- igraph2ggnetwork(graph, arrow.gap=0)

        # Plot graph
        p <- ggplot2::ggplot(n, ggplot2::aes_(x = ~x, y = ~y, xend = ~xend, yend = ~yend)) +
            ggnetwork::geom_edges(color = "grey75", alpha = 0.5, show.legend=FALSE) +
            ggnetwork::geom_nodes(ggplot2::aes_(size = ~Degree, color = ~Class)) +
            set_legend +
            ggplot2::scale_color_manual(values = palette) +
            add_nodelabel +
            ggnetwork::theme_blank()
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
#' @param ranked Logical indicating whether to treat third column of the edgelist (edge weights) as ranked values. Default: TRUE.
#'
#' @return A ggplot object containing the network.
#' @seealso
#'  \code{\link[networkD3]{igraph_to_networkD3}},\code{\link[networkD3]{forceNetwork}}
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
#' @examples
#' data(filt.se)
#' tfs <- sample(rownames(filt.se), size=50, replace=FALSE)
#' grn_edges <- grn_clr(filt.se, regulators = tfs)
#' p <- plot_grn(grn_edges, ranked=FALSE)
plot_grn <- function(edgelist_grn, show_labels = "tophubs", top_n_hubs = 5,
                     interactive = FALSE, layout = "kk", arrow.gap = 0.01,
                     ranked = TRUE) {
    requireNamespace("intergraph", quietly=TRUE)

    # Get degree of nodes and find hubs
    h <- get_hubs_grn(edgelist_grn, return_degree = TRUE, ranked=ranked)

    # Create a data frame of nodes and node attributes
    geneIDs <- unique(c(as.character(edgelist_grn[,1]), as.character(edgelist_grn[,2])))
    nod_at <- data.frame(Gene = geneIDs, stringsAsFactors = FALSE)
    nod_at$Class <- ifelse(nod_at$Gene %in% edgelist_grn[,1], "Regulator", "Target")
    nod_at$Degree <- h$Degree$Degree[h$Degree$Gene %in% nod_at$Gene]
    nod_at$isHub <- ifelse(nod_at$Gene %in% h$Hubs$Gene, TRUE, FALSE)
    nod_at <- nod_at[order(nod_at$Class, -nod_at$Degree), ]

    # Should the network be interactive?
    if(interactive) {
        graph <- igraph::graph_from_data_frame(d = edgelist_grn,
                                               vertices = nod_at, directed=TRUE)
        graph <- igraph::simplify(graph)
        graph_d3 <- networkD3::igraph_to_networkD3(graph, group = nod_at$Class)
        graph_d3$nodes <- merge(graph_d3$nodes, nod_at, by.x="name", by.y="Gene", sort = FALSE)
        my_color <- 'd3.scaleOrdinal() .domain(["Regulator", "Target"]) .range(["forestgreen", "orange"])'
        p <- networkD3::forceNetwork(Links = graph_d3$links, Nodes = graph_d3$nodes,
                                     Source = 'source', Target = 'target',
                                     NodeID = 'name', Group = 'group',
                                     colourScale = my_color,
                                     Nodesize = 'Degree', height=900, width=1200,
                                     opacity=1, zoom = TRUE, fontSize = 20, legend=TRUE)

    } else { #Static network
        # Define plotting parameters
        add_edges <- ggnetwork::geom_edges(color="grey60", alpha = 0.5,
                                           arrow = ggplot2::arrow(length = ggnetwork::unit(0.1, "lines"),
                                                                  type = "closed"),
                                           curvature = 0.1, show.legend = FALSE)

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
            tophubs <- nod_at[nod_at$isHub == TRUE, 1][seq_len(top_n_hubs)]
            nod_at$isTopHub <- ifelse(nod_at$Gene %in% tophubs, TRUE, FALSE)
            add_nodelabel <- ggnetwork::geom_nodelabel_repel(ggplot2::aes_(label = ~name),
                                                             color="azure4",
                                                             box.padding = ggnetwork::unit(1, "lines"),
                                                             data = function(x) { x[ x$isTopHub, ]},
                                                             show.legend = FALSE, max.overlaps=Inf)
        } else if(show_labels == "none") {
            add_nodelabel <- NULL
        } else {
            stop("Please, specify a valid option for 'show_labels'.")
        }
        # Create graph object
        graph <- igraph::graph_from_data_frame(d=edgelist_grn,
                                               vertices=nod_at, directed=TRUE)
        graph <- igraph::simplify(graph)
        n <- igraph2ggnetwork(graph, layout = layout, arrow.gap = arrow.gap)
        n$Class <- factor(n$Class, levels=c("Target", "Regulator"))

        # Plot graph
        p <- ggplot2::ggplot(n, ggplot2::aes_(x = ~x, y = ~y, xend = ~xend, yend = ~yend)) +
            add_edges +
            ggnewscale::new_scale_color() +
            ggnetwork::geom_nodes(ggplot2::aes_(color = ~Class, size = ~Degree, shape = ~Class)) +
            ggplot2::scale_color_manual(values = c("orange", "darkgreen")) +
            add_nodelabel +
            ggnetwork::theme_blank()
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
#' @examples
#' data(filt.se)
#' gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson",
#'                reportPDF = FALSE)
#' gcn_edges <- get_edge_list(gcn, module="brown", filter=TRUE,
#'                            method="min_cor")
#' hubs <- get_hubs_gcn(filt.se, gcn)
#' p <- plot_gcn(gcn_edges, gcn, hubs = hubs)
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
            palette <- custom_palette(1)[seq_len(nlevels(nod_at$Class))]
            scale_color <- ggplot2::scale_color_manual(values = palette)
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


