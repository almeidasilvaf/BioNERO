#' Calculate module preservation between two expression data sets using WGCNA's algorithm
#'
#' @param explist List of expression matrices with gene/probe names corresponding to column names and sample names corresponding to row names.
#' @param savePreservation Logical indicating whether to save module preservation into an R object or not. As the calculation of module preservation can take a long time, it is useful to save time in future analyses. Default is TRUE.
#' @param plot_all_stats Logical indicating whether to save all density and connectivity statistics in a PDF file or not. Default is FALSE.
#' @param calculateClusterCoeff Logical indicating if clustering coefficient statistics should be calculated. Only valid if \code{plot_all_stats} is TRUE. As it may take a long time, default is FALSE.
#'
#' @return A list object resulting from \code{WGCNA::modulePreservation} and .pdf files with plot statistics.
#' @seealso
#'  \code{\link[WGCNA]{modulePreservation}},\code{\link[WGCNA]{standardColors}}
#'  \code{\link[ggpubr]{ggscatter}},\code{\link[ggpubr]{ggarrange}},\code{\link[ggpubr]{ggexport}}
#'  \code{\link[ggplot2]{theme}},\code{\link[ggplot2]{margin}},\code{\link[ggplot2]{geom_abline}}
#' @rdname modPres_WGCNA
#' @export
#' @importFrom WGCNA modulePreservation standardColors
#' @importFrom ggpubr ggscatter ggarrange ggexport
#' @importFrom ggplot2 theme element_text geom_hline
modPres_WGCNA <- function(explist, ref_net,
                          savePreservation = TRUE, plot_all_stats = FALSE,
                          calculateClusterCoeff = FALSE) {

    # Keep only genes/probes present in both data sets
    overlap <- as.character(intersect(colnames(explist[[1]]), colnames(explist[[2]])))
    explist[[1]] <- explist[[1]][, overlap]
    explist[[2]] <- explist[[2]][, overlap]

    # Set parameters for network reconstruction
    net_type <- ref_net$params$net_type
    cor_method <- ref_net$params$cor_method

    # Create multiExpr object
    multiExpr <- list(ref = list(data = explist[[1]]),
                      test = list(data = explist[[2]]))

    # Create vector of module assignments
    moduleColors <- ref_net$genes_and_modules$Modules
    names(moduleColors) <- ref_net$genes_and_modules$Genes

    multiColor <- list(ref = moduleColors)

    # Calculate module preservation
    if(cor_method == "pearson") {
        pres <- WGCNA::modulePreservation(multiExpr, multiColor,
                                          referenceNetworks = 1,
                                          nPermutations = 200,
                                          randomSeed = 1,
                                          quickCor = 0,
                                          verbose = 3, networkType=net_type,
                                          calculateClusterCoeff = calculateClusterCoeff)
    } else if(cor_method == "spearman") {
        pres <- WGCNA::modulePreservation(multiExpr, multiColor,
                                          referenceNetworks = 1,
                                          nPermutations = 200,
                                          randomSeed = 1,
                                          quickCor = 0,
                                          corOptions = list(use="p", method="spearman"),
                                          verbose = 3, networkType=net_type,
                                          calculateClusterCoeff = calculateClusterCoeff)
    } else if(cor_method == "biweight") {
        pres <- WGCNA::modulePreservation(multiExpr, multiColor,
                                          referenceNetworks = 1,
                                          nPermutations = 200,
                                          randomSeed = 1,
                                          quickCor = 0, corFnc = "bicor",
                                          verbose = 3, networkType=net_type,
                                          calculateClusterCoeff = calculateClusterCoeff)
    } else {
        stop("Please, specify a valid correlation method.")
    }

    if(savePreservation == TRUE) {
        save(pres, file = "modulePreservation.RData");
    }

    # Isolate the observed statistics and their Z-scores
    ref <- 1
    test <- 2
    statsObs <- cbind(pres$quality$observed[[ref]][[test]][, -1], pres$preservation$observed[[ref]][[test]][, -1])
    statsZ <- cbind(pres$quality$Z[[ref]][[test]][, -1], pres$preservation$Z[[ref]][[test]][, -1])


    # Module labels and module sizes are also contained in the results
    modColors <- rownames(pres$preservation$observed[[ref]][[test]])
    moduleSizes <- pres$preservation$Z[[ref]][[test]][, 1]

    # Do not consider grey and gold modules
    plotMods <- !(modColors %in% c("grey", "gold"))

    # Text labels for points
    text <- modColors[plotMods]

    # Auxiliary convenience variable
    plotData <- cbind(pres$preservation$observed[[ref]][[test]][, 2], pres$preservation$Z[[ref]][[test]][, 2])

    # Main titles for the plot
    mains <- c("Preservation Median rank", "Preservation Zsummary");

    # Plot preservation median rank
    dplot1 <- data.frame(x=moduleSizes[plotMods], y=plotData[plotMods, 1], label=text, stringsAsFactors=FALSE)
    p1 <- ggpubr::ggscatter(dplot1, x = "x", y = "y", size=4,
                            color="black", fill=text, shape=21,
                            label = "label", repel = TRUE,
                            xlab = "Module size", ylab = mains[1],
                            title = mains[1], font.title = c(17, "bold")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    # Plot preservation Z-summary
    dplot2 <- data.frame(x=moduleSizes[plotMods], y=plotData[plotMods, 2], label=text, stringsAsFactors=FALSE)
    p2 <- ggscatter(dplot2, x = "x", y = "y", size=4,
                    color="black", fill=text, shape=21,
                    label = "label", repel = TRUE,
                    xlab = "Module size", ylab = mains[2],
                    title = mains[2], font.title = c(17, "bold")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::geom_hline(yintercept = 2, colour = "blue", linetype = 2) +
        ggplot2::geom_hline(yintercept = 10, colour = "forestgreen", linetype = 2)

    fig1 <- ggpubr::ggarrange(p1, p2, ncol=2, nrow=1)
    ggpubr::ggexport(fig1, filename = "Module_preservation.pdf", height=6, width=12)


    if(plot_all_stats == TRUE) {

        # Re-initialize module color labels and sizes
        modColors <- rownames(statsZ)
        moduleSizes <- pres$quality$Z[[ref]][[test]][, 1];

        # Exclude improper modules
        plotMods <- !(modColors %in% c("grey", "gold"));

        # Create numeric labels for each module
        labs <- match(modColors[plotMods], WGCNA::standardColors(50));

        # Create a list of plots per column of statsZ
        df <- cbind(moduleSizes[plotMods], statsZ[plotMods, ])
        colnames(df)[1] <- "module_size"
        ylabs <- colnames(df)[2:ncol(df)]
        df$labels <- rownames(df)

        plots <- lapply(ylabs, function(x) ggpubr::ggscatter(df, x = "module_size", y = x, size=4,
                                                             color="black", fill=labels, shape=21,
                                                             label = "labels", repel = TRUE,
                                                             xlab = "Module size", ylab = x,
                                                             title = x, font.title = c(17, "bold")) +
                            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
                            ggplot2::geom_hline(yintercept = 0) +
                            ggplot2::geom_hline(yintercept = 2, colour = "blue", linetype = 2) +
                            ggplot2::geom_hline(yintercept = 10, colour = "forestgreen", linetype = 2))

        if(calculateClusterCoeff == FALSE) {
            pl1 <- plots[1:4]
            pl2 <- plots[5:8]
            pl3 <- plots[9:12]
            pl4 <- plots[13:16]
            pl5 <- plots[17:19]
            mpage1 <- ggpubr::ggarrange(plotlist = pl1)
            mpage2 <- ggpubr::ggarrange(plotlist = pl2)
            mpage3 <- ggpubr::ggarrange(plotlist = pl3)
            mpage4 <- ggpubr::ggarrange(plotlist = pl4)
            mpage5 <- ggpubr::ggarrange(plotlist = pl5)
            final_list <- list(mpage1, mpage2, mpage3, mpage4, mpage5)

            multipage <- ggpubr::ggarrange(plotlist = final_list, ncol=1, nrow=1)
            ggpubr::ggexport(multipage, filename="all_Zstats.pdf", height = 9, width=9)
        } else {
            pl1 <- plots[1:4]
            pl2 <- plots[5:8]
            pl3 <- plots[9:12]
            pl4 <- plots[13:16]
            pl5 <- plots[17:20]
            pl6 <- plots[21:22]
            mpage1 <- ggpubr::ggarrange(plotlist = pl1)
            mpage2 <- ggpubr::ggarrange(plotlist = pl2)
            mpage3 <- ggpubr::ggarrange(plotlist = pl3)
            mpage4 <- ggpubr::ggarrange(plotlist = pl4)
            mpage5 <- ggpubr::ggarrange(plotlist = pl5)
            mpage6 <- ggpubr::ggarrange(plotlist = pl6)
            final_list <- list(mpage1, mpage2, mpage3, mpage4, mpage5, mpage6)

            multipage <- ggpubr::ggarrange(plotlist = final_list, ncol=1, nrow=1)
            ggpubr::ggexport(multipage, filename = "all_Zstats.pdf", height = 9, width=9)
        }
    }
    return(pres)
}


#' Calculate module preservation between two expression data sets using NetRep's algorithm
#'
#' @param explist List of expression matrices with gene/probe names corresponding to column names and sample names corresponding to row names.
#' @param ref_net Reference network object returned by the function \code{exp2net}.
#' @param test_net Test network object returned by the function \code{exp2net}.
#' @param nPerm Number of permutations. Default: 1000
#' @param nThreads Number of threads to be used for parallel computing. Default: 2
#' @return Output list from \code{NetRep::modulePreservation} and a message in user's standard output stating which modules are preserved.
#' @seealso
#'  \code{\link[NetRep]{modulePreservation}}
#' @rdname modPres_netrep
#' @export
#' @importFrom NetRep modulePreservation
modPres_netrep <- function(explist, ref_net = NULL, test_net = NULL,
                           nPerm = 1000, nThreads=2) {

    # Keep only genes/probes present in both data sets
    overlap <- as.character(intersect(colnames(explist[[1]]), colnames(explist[[2]])))
    explist[[1]] <- explist[[1]][, overlap]
    explist[[2]] <- explist[[2]][, overlap]

    # Set data set names
    if(is.null(names(explist))) {
        data_names <- c("cohort1", "cohort2")
    } else {
        data_names <- names(explist)
    }

    # Create correlation list
    correlation_list <- list(as.matrix(ref_net$correlation_matrix),
                             as.matrix(test_net$correlation_matrix))
    names(correlation_list) <- data_names

    # Create network list
    network_list <- list(as.matrix(ref_net$adjacency_matrix),
                         as.matrix(test_net$adjacency_matrix))
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
        network=network_list, data=explist, correlation=correlation_list,
        moduleAssignments=modAssignments,
        discovery=data_names[1], test=data_names[2],
        nPerm=nPerm, nThreads=nThreads
    )

    # Get preserved modules (p < 0.05 for all statistics)
    max_pval <- apply(pres2$p.value, 1, max)
    preservedmodules <- names(max_pval[max_pval < 0.05])
    cat("We found", length(preservedmodules), "preserved modules:", preservedmodules)

    return(pres)
}



#' Calculate network preservation between two expression data sets
#'
#' @param explist List of expression matrices with gene/probe names corresponding to column names and sample names corresponding to row names.
#' @param ref_net Reference network object returned by the function \code{exp2net}.
#' @param test_net Test network object returned by the function \code{exp2net}.
#' @param algorithm Module preservation algorithm to be used. One of 'netrep' (default, permutation-based) or WGCNA.
#' @param diffIDs Logical indicating whether the different data sets have different gene IDs (e.g., for different species). Default: FALSE
#' @param correspondence Data frame containing gene IDs for the reference data set in the column 1 and gene IDs for the test data set in column 2. Only required if diffIDs is TRUE.
#' @param nPerm Number of permutations. Default: 1000
#' @param nThreads Number of threads to be used for parallel computing. Default: 2
#' @param savePreservation Logical indicating whether to save module preservation into an R object or not. As the calculation of module preservation can take a long time, it is useful to save time in future analyses. Default is TRUE.
#' @param plot_all_stats Logical indicating whether to save all density and connectivity statistics in a PDF file or not. Default is FALSE.
#' @param calculateClusterCoeff Logical indicating if clustering coefficient statistics should be calculated. Only valid if \code{plot_all_stats} is TRUE. As it may take a long time, default is FALSE.
#' @return A list containing the preservation statistics. The list will be different depending on the chosen algorithm. See \code{WGCNA::modulePreservation} or \code{NetRep::modulePreservation} for more info.
#' @rdname module_preservation
#' @export
module_preservation <- function(explist, ref_net = NULL, test_net = NULL,
                                algorithm="netrep",
                                diffIDs=FALSE,
                                correspondence=NULL,
                                nPerm = 1000, nThreads=2,
                                savePreservation = TRUE,
                                plot_all_stats = FALSE,
                                calculateClusterCoeff = FALSE) {

    if(diffIDs == TRUE) {
        if(is.null(correspondence)) {
            stop("Please, input a data frame of correspondence between IDs from the different data sets.")
        }

        # Set IDs of test data set to be equal to IDs of reference data set
        colnames(explist[[2]]) <- colnames(explist[[1]])
        colnames(explist[[2]]) <- correspondence[correspondence[,2] %in%
                                                     colnames(explist[[2]]), 1]
    }

    if(algorithm == "netrep") {
        pres <- modPres_netrep(explist = explist,
                               ref_net = ref_net,
                               test_net = test_net,
                               nPerm = nPerm, nThreads = nThreads)

    } else if(algorithm == "WGCNA") {
        pres <- modPres_WGCNA(explist = explist,
                              ref_net = ref_net,
                              savePreservation = savePreservation,
                              plot_all_stats = plot_all_stats,
                              calculateClusterCoeff = calculateClusterCoeff)
    } else {
        stop("Please, specify a valid algorithm. One of 'netrep' or 'WGCNA'.")
    }

    return(pres)
}


