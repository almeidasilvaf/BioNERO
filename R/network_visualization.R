#' Detect communities in a network
#'
#' @param edgelist Data frame containing the network as an edge list.
#' First column must be node 1 and second column must be node 2.
#' Additional columns will be interpreted as edge attributes and will
#' be modified by this function.
#' @param method igraph function to be used for community detection.
#' Available functions are cluster_infomap, cluster_edge_betweenness,
#' cluster_fast_greedy, cluster_walktrap, cluster_spinglass,
#' cluster_leading_eigen, cluster_louvain, and cluster_label_prop.
#' Default is cluster_infomap.
#' @param directed Logical indicating whether the network is directed (GRN only)
#' or not (GCN and PPI networks). Default: TRUE.
#' @return A data frame containing node names in the first column, and
#' communities to which nodes belong in the second column.
#' @seealso
#'  \code{\link[igraph]{cluster_infomap}},
#'  \code{\link[igraph]{cluster_edge_betweenness}},
#'  \code{\link[igraph]{cluster_fast_greedy}},
#'  \code{\link[igraph]{cluster_walktrap}},
#'  \code{\link[igraph]{cluster_spinglass}},
#'  \code{\link[igraph]{cluster_leading_eigen}},
#'  \code{\link[igraph]{cluster_louvain}},
#'  \code{\link[igraph]{cluster_label_prop}}
#' @rdname detect_communities
#' @author Fabricio Almeida-Silva
#' @export
#' @importFrom igraph graph.data.frame simplify cluster_infomap cluster_edge_betweenness cluster_fast_greedy cluster_walktrap cluster_spinglass cluster_leading_eigen cluster_louvain cluster_label_prop
#' @examples
#' data(filt.se)
#' tfs <- sample(rownames(filt.se), size=50, replace=FALSE)
#' grn_edges <- grn_infer(filt.se, method = "clr", regulators = tfs)
#' com <- detect_communities(grn_edges, directed=TRUE)
detect_communities <- function(edgelist,
                               method = igraph::cluster_infomap,
                               directed = TRUE) {
    graph <- igraph::graph.data.frame(edgelist, directed = directed)
    graph <- igraph::simplify(graph)
    com <- method(graph, modularity = FALSE)
    df_com <- as.data.frame(list(names = com$names, mem = com$membership))
    return(df_com)
}


#' Plot protein-protein interaction network from edge list
#'
#' @param edgelist_int Data frame containing the edge list for the PPI network.
#' First column is the protein 1 and second column is the protein 2.
#' All other columns are interpreted as edge attributes.
#' @param color_by How should nodes be colored? It must be either "community" or
#' a 2-column data frame containing proteins in the first column and a
#' custom annotation in the second column. If "community", a clustering
#' algorithm will be applied. Default: "community".
#' @param clustering_method igraph function to be used for community detection.
#' Available functions are cluster_infomap, cluster_edge_betweenness,
#' cluster_fast_greedy, cluster_walktrap, cluster_spinglass,
#' cluster_leading_eigen, cluster_louvain, and cluster_label_prop.
#' Default is cluster_infomap.
#' @param show_labels Character indicating which nodes will be labeled.
#' One of "all", "allhubs", "tophubs", or "none".
#' @param top_n_hubs Number of top hubs to be labeled. It is only valid
#' if \code{show_labels} equals "tophubs". Default is 5.
#' @param interactive Logical indicating whether the network should be
#' interactive or not. Default is FALSE.
#' @param add_color_legend Logical indicating whether to add a color legend
#' for nodes. Default: TRUE.
#' @param dim_interactive Numeric vector with width and height of window
#' for interactive plotting. Default: c(600,600).
#'
#' @return A ggplot object.
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
#' @import intergraph
#' @examples
#' ppi_edges <- igraph::get.edgelist(igraph::barabasi.game(n=50, directed=FALSE))
#' p <- plot_ppi(ppi_edges, add_color_legend = FALSE)
plot_ppi <- function(edgelist_int, color_by = "community",
                     clustering_method = igraph::cluster_infomap,
                     show_labels = "tophubs",
                     top_n_hubs = 5, interactive = FALSE,
                     add_color_legend=TRUE, dim_interactive = c(600,600)) {
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
        d <- dim_interactive
        p <- networkD3::forceNetwork(Links = graph_d3$links, Nodes = graph_d3$nodes,
                                     Source = 'source', Target = 'target',
                                     NodeID = 'name', Group = 'group',
                                     Nodesize = 'Degree', height=d[2], width=d[1],
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
            tophubs <- nod_at[nod_at$isHub == TRUE, 1][seq_len(top_n_hubs)]
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
        n <- ggnetwork::ggnetwork(graph, arrow.gap = 0)

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
#' @param edgelist_grn Data frame containing the edge list for the GRN network.
#' First column is the TF and second column is the target gene.
#' All other columns are interpreted as edge attributes.
#' @param show_labels Character indicating which nodes will be labeled.
#' One of "all", "allhubs", "tophubs", or "none".
#' @param top_n_hubs Number of top hubs to be labeled. It is only valid
#' if \code{show_labels} equals "tophubs". Default is 5.
#' @param interactive Logical indicating whether the network should be
#' interactive or not. Default is FALSE.
#' @param layout igraph function for the network layout. One of
#' with_dh, with_drl, with_gem, with_lgl, with_fr, with_graphopt,
#' with_kk and with_mds. Default is with_kk.
#' @param arrow.gap Numeric indicating the distance between nodes and arrows.
#' Default is 0.2.
#' @param ranked Logical indicating whether to treat third column of
#' the edge list (edge weights) as ranked values. Default: TRUE.
#' @param dim_interactive Numeric vector with width and height of window
#' for interactive plotting. Default: c(600,600).
#'
#' @return A ggplot object containing the network.
#' @seealso
#'  \code{\link[networkD3]{igraph_to_networkD3}},
#'  \code{\link[networkD3]{forceNetwork}}
#'  \code{\link[ggnetwork]{geom_edges}}, \code{\link[ggnetwork]{geom_nodes}},
#'  \code{\link[ggnetwork]{geom_nodetext}},\code{\link[ggnetwork]{theme_blank}},
#'  \code{\link[ggnetwork]{geom_nodetext_repel}}
#'  \code{\link[ggnewscale]{new_scale}}
#' @rdname plot_grn
#' @author Fabricio Almeida-Silva
#' @export
#' @import intergraph
#' @importFrom igraph graph_from_data_frame degree vcount
#' @importFrom networkD3 igraph_to_networkD3 forceNetwork
#' @importFrom ggplot2 ggplot aes_ arrow scale_color_manual guides
#' @importFrom ggnetwork geom_edges unit geom_nodes geom_nodetext theme_blank geom_nodelabel_repel
#' @importFrom ggnewscale new_scale_color
#' @examples
#' data(filt.se)
#' tfs <- sample(rownames(filt.se), size=50, replace=FALSE)
#' grn_edges <- grn_infer(filt.se, method ="clr", regulators = tfs)
#' p <- plot_grn(grn_edges, ranked=FALSE)
plot_grn <- function(edgelist_grn, show_labels = "tophubs", top_n_hubs = 5,
                     interactive = FALSE,
                     layout = igraph::with_kk, arrow.gap = 0.01,
                     ranked = TRUE, dim_interactive = c(600,600)) {
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
        d <- dim_interactive
        p <- networkD3::forceNetwork(Links = graph_d3$links, Nodes = graph_d3$nodes,
                                     Source = 'source', Target = 'target',
                                     NodeID = 'name', Group = 'group',
                                     colourScale = my_color,
                                     Nodesize = 'Degree', height=d[2], width=d[1],
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
        n <- ggnetwork::ggnetwork(graph, layout = layout(), arrow.gap = arrow.gap)
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
#' @param edgelist_gcn Data frame containing the edge list for the GCN.
#' The edge list can be generated with \code{get_edge_list()}.
#' @param net List object returned by \code{exp2net}.
#' @param color_by How should nodes be colored? It must be either "module"
#' (nodes will have the colors of their modules) or a 2-column data frame
#' containing genes in the first column and a custom gene annotation
#' in the second column. Default: "module".
#' @param hubs Data frame containing hub genes in the first column,
#' their modules in the second column, and intramodular connectivity in
#' the third column.
#' @param show_labels Character indicating which nodes will be labeled.
#' One of "all", "allhubs", "tophubs", or "none". Default: tophubs.
#' @param top_n_hubs Number of top hubs to be labeled. It is only valid
#' if \code{show_labels} equals "tophubs". Default is 5.
#' @param interactive Logical indicating whether the network should be
#' interactive or not. Default is FALSE.
#' @param dim_interactive Numeric vector with width and height of window
#' for interactive plotting. Default: c(600,600).
#'
#' @return A ggplot object.
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
#' @import intergraph
#' @examples
#' data(filt.se)
#' gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson")
#' gcn_edges <- get_edge_list(gcn, module="brown", filter=TRUE,
#'                            method="min_cor")
#' hubs <- get_hubs_gcn(filt.se, gcn)
#' p <- plot_gcn(gcn_edges, gcn, hubs = hubs)
plot_gcn <- function(edgelist_gcn, net, color_by="module", hubs = NULL,
                     show_labels = "tophubs", top_n_hubs = 5,
                     interactive = FALSE, dim_interactive = c(600,600)) {
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
    nod_at <- merge(nod_at, gene_annotation, by = 1)
    names(nod_at)[2] <- "Class"
    nod_at$Class <- as.factor(nod_at$Class)
    nod_at$Degree <- kIN$kWithin[rownames(kIN) %in% nod_at$Gene]
    nod_at$isHub <- ifelse(nod_at$Gene %in% hubs[,1], TRUE, FALSE)
    nod_at <- nod_at[order(nod_at$Class, -nod_at$Degree), ]

    # Should the network be interactive?
    if(interactive) {
        graph <- igraph::simplify(igraph::graph_from_data_frame(d = edgelist_gcn, vertices = nod_at, directed=FALSE))
        graph_d3 <- networkD3::igraph_to_networkD3(graph, group = nod_at$Class)
        graph_d3$nodes <- merge(graph_d3$nodes, nod_at, by.x="name", by.y="Gene", sort = FALSE)
        d <- dim_interactive
        p <- networkD3::forceNetwork(Links = graph_d3$links, Nodes = graph_d3$nodes,
                                     Source = 'source', Target = 'target',
                                     NodeID = 'name', Group = 'group',
                                     Nodesize = 'Degree', height=d[2], width=d[1],
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
            tophubs <- nod_at[nod_at$isHub == TRUE, 1][seq_len(top_n_hubs)]
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
        n <- ggnetwork::ggnetwork(graph, arrow.gap=0)

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


