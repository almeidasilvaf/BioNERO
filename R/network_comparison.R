
#' Parse orthogroups identified by OrthoFinder
#'
#' This function converts the orthogroups file named \strong{Orthogroups.tsv} to
#' a 3-column data frame that can be interpreted by BioNERO.
#'
#' @param file_path Path to Orthogroups/Orthogroups.tsv file generated
#' by OrthoFinder.
#'
#' @author Fabricio Almeida-Silva
#' @return A 3-column data frame with orthogroups, species IDs and
#' gene IDs, respectively.
#' @importFrom reshape2 melt
#' @export
#' @rdname parse_orthofinder
#' @examples
#' path <- system.file("extdata", "Orthogroups.tsv", package = "BioNERO")
#' og <- parse_orthofinder(path)
parse_orthofinder <- function(file_path = NULL) {
  of <- read.csv(file_path, sep="\t")
  melt_of <- reshape2::melt(of, id.vars="Orthogroup")
  s <- strsplit(melt_of$value, split = ", ")
  final_of <- data.frame(
      Orthogroup = rep(melt_of$Orthogroup, vapply(s, length, numeric(1))),
      Species = rep(melt_of$variable, vapply(s, length, numeric(1))),
      Gene = unlist(s)
  )
  return(final_of)
}


#' Collapse gene-level expression data to orthogroup level
#'
#' For a given list of expression data, this function replaces genes with
#' their corresponding orthogroups to allow inter-species comparisons.
#'
#' @param explist List of expression data frames or SummarizedExperiment objects.
#' @param og Data frame of 3 columns corresponding to orthogroup,
#' species ID, and gene ID, respectively.
#' Species IDs must be the same as the names of the expression list.
#' @param summarize Centrality measure to summarize multiple paralogous genes
#' in the same orthogroup.
#' One of "median" or "mean". Default: "median".
#' @return List of expression data frames for each species with expression
#' summarized at the orthogroup level.
#'
#' @export
#' @importFrom stats aggregate
#' @rdname exp_genes2orthogroups
#' @examples
#' \donttest{
#' data(og.zma.osa)
#' data(zma.se)
#' data(osa.se)
#' explist <- list(zma = zma.se,
#'                 osa = osa.se)
#' og <- og.zma.osa
#' exp_ortho <- exp_genes2orthogroups(explist, og, summarize = "mean")
#' }
exp_genes2orthogroups <- function(explist = NULL, og = NULL,
                                  summarize = "median") {

    colnames(og) <- c("Family", "Species", "Gene")
    exp <- handleSElist(explist)
    exp <- exp[order(names(exp))]
    og_list <- split(og, og$Species)
    og_list <- og_list[order(names(og_list))]
    og_exp <- lapply(seq_along(exp), function(x) {
        set <- merge(exp[[x]], og_list[[x]], by.x = "row.names", by.y = 3)
        set <- set[, -c(1, ncol(set))]
        set <- suppressWarnings(
            aggregate(set, by = list(set$Family), FUN = summarize)
        )
        rownames(set) <- set[,1]
        set <- set[, -c(1, ncol(set))]
        return(set)
    })
    names(og_exp) <- names(exp)

    return(og_exp)
}

#' Calculate module preservation between two expression data sets using WGCNA's algorithm
#'
#' @param explist List of expression data frames or SummarizedExperiment objects.
#' @param ref_net Reference network object returned by the function \code{exp2net}.
#' @param nPerm Number of permutations for the module preservation statistics.
#' It must be greater than 1. Default: 200.
#'
#' @return A ggplot object with module preservation statistics.
#' @rdname modPres_WGCNA
#' @export
#' @importFrom WGCNA standardColors
#' @importFrom ggplot2 theme geom_hline geom_point
#' @importFrom ggrepel geom_text_repel
#' @importFrom patchwork wrap_plots
#' @examples
#' \donttest{
#' set.seed(1)
#' data(og.zma.osa)
#' data(zma.se)
#' data(osa.se)
#' explist <- list(Zma = zma.se, Osa = osa.se)
#' og <- og.zma.osa
#' exp_ortho <- exp_genes2orthogroups(explist, og, summarize = "mean")
#' exp_ortho <- lapply(exp_ortho, function(x) filter_by_variance(x, n=1500))
#' # Previously calculated power
#' powers <- c(13, 15)
#' gcn_osa <- exp2gcn(exp_ortho$Osa, net_type = "signed hybrid",
#'                    SFTpower = powers[1], cor_method = "pearson")
#' explist <- exp_ortho
#' ref_net <- gcn_osa
#' # 5 permutations for demonstration purposes
#' pres_wgcna <- modPres_WGCNA(explist, ref_net, nPerm=5)
#' }
modPres_WGCNA <- function(explist, ref_net, nPerm = 200) {

    explist <- handleSElist(explist)
    explist <- lapply(explist, function(x) return(t(x)))

    # Set parameters for network reconstruction
    net_type <- ref_net$params$net_type
    cor_method <- ref_net$params$cor_method

    # Create multiExpr object
    multiExpr <- list(
        ref = list(data = explist[[1]]), test = list(data = explist[[2]])
    )

    # Create vector of module assignments
    moduleColors <- ref_net$genes_and_modules$Modules
    names(moduleColors) <- ref_net$genes_and_modules$Genes

    multiColor <- list(ref = moduleColors)

    # Calculate module preservation
    if(cor_method == "pearson") {
      corOptions <- "use = 'p'"
      corFnc <- "cor"
    } else if(cor_method == "spearman") {
      corOptions <- list(use = 'p', method="spearman")
      corFnc <- "cor"
    } else if(cor_method == "biweight") {
      corOptions <- list(use = 'p', maxPOutliers = 0.05)
      corFnc <- "bicor"
    } else {
      stop("Please, specify a valid correlation method.")
    }

    pres <- WGCNA::modulePreservation(
        multiExpr, multiColor, referenceNetworks = 1, nPermutations = nPerm,
        randomSeed = 1, quickCor = 0, corFnc = corFnc, corOptions = corOptions,
        verbose = 0, networkType = net_type, savePermutedStatistics = FALSE,
        plotInterpolation = FALSE
    )

    # Isolate the observed statistics and their Z-scores
    ref <- 1
    test <- 2
    statsObs <- cbind(pres$quality$observed[[ref]][[test]][, -1],
                      pres$preservation$observed[[ref]][[test]][, -1])
    statsZ <- cbind(pres$quality$Z[[ref]][[test]][, -1],
                    pres$preservation$Z[[ref]][[test]][, -1])

    # Module labels and module sizes are also contained in the results
    modColors <- rownames(pres$preservation$observed[[ref]][[test]])
    moduleSizes <- pres$preservation$Z[[ref]][[test]][, 1]

    # Do not consider grey and gold modules
    plotMods <- !(modColors %in% c("grey", "gold"))

    # Text labels for points
    text <- modColors[plotMods]

    # Auxiliary convenience variable
    plotData <- cbind(pres$preservation$observed[[ref]][[test]][, 2],
                      pres$preservation$Z[[ref]][[test]][, 2])


    # Plot preservation median rank
    dplot1 <- data.frame(
        x = moduleSizes[plotMods],
        y = plotData[plotMods, 1],
        label = as.factor(text)
    )

    p1 <- ggplot(dplot1, aes(x = .data$x, y = .data$y)) +
        geom_point(fill = text, color = "black", shape = 21) +
        geom_text_repel(aes(label = .data$label)) +
        labs(
            x = "Module size", y = "Median rank",
            title = "Preservation median rank"
        ) +
        theme_bw()

    # Plot preservation Z-summary
    dplot2 <- data.frame(
        x = moduleSizes[plotMods],
        y = plotData[plotMods, 2],
        label = text
    )

    p2 <- ggplot(dplot2, aes(x = .data$x, y = .data$y)) +
        geom_point(fill = text, color = "black", shape = 21) +
        geom_text_repel(aes(label = .data$label)) +
        labs(
            x = "Module size", y = expression("Z"[summary]),
            title = expression("Preservation Z"[summary])
        ) +
        theme_bw() +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_hline(yintercept = 2, colour = "blue", linetype = 2) +
        geom_hline(yintercept = 10, colour = "forestgreen", linetype = 2)

    fig1 <- patchwork::wrap_plots(p1, p2, ncol = 2)
    return(fig1)
}


#' Calculate module preservation between two expression data sets using NetRep's algorithm
#'
#' @param explist List of expression data frames or SummarizedExperiment objects.
#' @param ref_net Reference network object returned by the function \code{exp2net}.
#' @param test_net Test network object returned by the function \code{exp2net}.
#' @param nPerm Number of permutations. Default: 1000
#' @param nThreads Number of threads to be used for parallel computing.
#' Default: 1
#' @return Output list from \code{NetRep::modulePreservation} and a message in
#' user's standard output stating which modules are preserved.
#' @seealso
#'  \code{\link[NetRep]{modulePreservation}}
#' @rdname modPres_netrep
#' @export
#' @importFrom NetRep modulePreservation
#' @examples
#' \donttest{
#' set.seed(1)
#' data(og.zma.osa)
#' data(zma.se)
#' data(osa.se)
#' og <- og.zma.osa
#' exp_ortho <- exp_genes2orthogroups(explist, og, summarize = "mean")
#' exp_ortho <- lapply(exp_ortho, function(x) filter_by_variance(x, n=1500))
#' # Previously calculated SFT powers
#' powers <- c(13, 15)
#' gcn_osa <- exp2gcn(exp_ortho$osa, net_type = "signed hybrid",
#'                    SFTpower = powers[1], cor_method = "pearson")
#' gcn_zma <- exp2gcn(exp_ortho$zma, net_type = "signed hybrid",
#'                    SFTpower = powers[2], cor_method = "pearson")
#' explist <- exp_ortho
#' ref_net <- gcn_osa
#' test_net <- gcn_zma
#' # 10 permutations for demonstration purposes
#' pres_netrep <- modPres_netrep(explist, ref_net, test_net,
#'                               nPerm=10, nThreads = 2)
#' }
#'
modPres_netrep <- function(explist, ref_net = NULL, test_net = NULL,
                           nPerm = 1000, nThreads = 1) {
    explist <- handleSElist(explist)
    explist <- lapply(explist, function(x) return(t(x)))

    # Set data set names
    if(is.null(names(explist))) {
        data_names <- c("cohort1", "cohort2")
    } else {
        data_names <- names(explist)
    }

    # Create correlation list
    correlation_list <- list(
        as.matrix(ref_net$correlation_matrix),
        as.matrix(test_net$correlation_matrix)
    )
    names(correlation_list) <- data_names

    # Create network list
    network_list <- list(
        as.matrix(ref_net$adjacency_matrix),
        as.matrix(test_net$adjacency_matrix)
    )
    names(network_list) <- data_names

    # Set background label
    if("grey" %in% ref_net$moduleColors) {
        backgroundLabel <- "grey"
    } else {
        backgroundLabel <- 0
    }

    # Set module assignments
    modAssignments <- ref_net$genes_and_modules$Modules
    names(modAssignments) <- ref_net$genes_and_modules$Genes

    # Calculate preservation statistics
    pres <- NetRep::modulePreservation(
        network = network_list, data = explist, correlation = correlation_list,
        moduleAssignments = modAssignments,
        discovery = data_names[1], test = data_names[2],
        nPerm = nPerm, nThreads = nThreads
    )

    # Get preserved modules (p < 0.05 for all statistics)
    max_pval <- apply(pres$p.value, 1, max)
    preservedmodules <- names(max_pval[max_pval < 0.05])
    if(length(preservedmodules) > 0) {
        message(length(preservedmodules), " modules in ", data_names[1],
                " were preserved in ", data_names[2], ":", "\n",
                toString(preservedmodules))
    } else {
        message("None of the modules in ", data_names[1],
                " were preserved in ", data_names[2], ".")
    }
    return(pres)
}



#' Calculate network preservation between two expression data sets
#'
#' @param explist List of SummarizedExperiment objects or expression data frames
#' with genes (or orthogroups) in row names and samples in column names.
#' @param ref_net Reference network object returned by
#' the function \code{exp2gcn}.
#' @param test_net Test network object returned by the function \code{exp2gcn}.
#' @param algorithm Module preservation algorithm to be used. One of 'netrep'
#' (default, permutation-based) or WGCNA.
#' @param nPerm Number of permutations. Default: 1000
#' @param nThreads Number of threads to be used for parallel computing.
#' Default: 1
#'
#' @return A list containing the preservation statistics (netrep) or a ggplot
#' object with preservation statistics.
#' See \code{WGCNA::modulePreservation} or \code{NetRep::modulePreservation}
#' for more info.
#' @rdname module_preservation
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' data(og.zma.osa)
#' data(zma.se)
#' data(osa.se)
#' og <- og.zma.osa
#' exp_ortho <- exp_genes2orthogroups(explist, og, summarize = "mean")
#' exp_ortho <- lapply(exp_ortho, function(x) filter_by_variance(x, n=1500))
#' # Previously calculated SFT powers
#' powers <- c(13, 15)
#' gcn_osa <- exp2gcn(exp_ortho$osa, net_type = "signed hybrid",
#'                    SFTpower = powers[1], cor_method = "pearson")
#' gcn_zma <- exp2gcn(exp_ortho$zma, net_type = "signed hybrid",
#'                    SFTpower = powers[2], cor_method = "pearson")
#' explist <- exp_ortho
#' ref_net <- gcn_osa
#' test_net <- gcn_zma
#' # 10 permutations for demonstration purposes
#' pres <- module_preservation(explist, ref_net, test_net, nPerm=10)
#' }
#'
module_preservation <- function(explist, ref_net = NULL, test_net = NULL,
                                algorithm = "netrep",
                                nPerm = 1000, nThreads = 1) {

    if(algorithm == "netrep") {
        pres <- modPres_netrep(
            explist = explist,
            ref_net = ref_net,
            test_net = test_net,
            nPerm = nPerm, nThreads = nThreads
        )

    } else if(algorithm == "WGCNA") {
        pres <- modPres_WGCNA(
            explist = explist,
            ref_net = ref_net,
            nPerm = nPerm
        )
    } else {
        stop("Please, specify a valid algorithm. One of 'netrep' or 'WGCNA'.")
    }

    return(pres)
}


