#' Calculate module preservation between two expression data sets
#'
#' @param exp_ref Data frame of reference gene expression data set, with gene IDs in row names and sample names in column names.
#' @param exp_test Data frame of test gene expression data set, with gene IDs in row names and sample names in column names.
#' @param moduleColors Vector of module assignment returned by \code{exp2net}.
#' @param savePreservation Logical indicating whether to save module preservation into an R object or not. As the calculation of module preservation can take a long time, it is useful to save time in future analyses. Default is TRUE.
#' @param plot_all_stats Logical indicating whether to save all density and connectivity statistics in a PDF file or not. Default is FALSE.
#'
#' @return Module preservation statistics.
#' @rdname module_preservation
#' @seealso
#'  \code{\link[WGCNA]{modulePreservation}}
#' @export
#' @importFrom WGCNA modulePreservation
module_preservation <- function(exp_ref, exp_test, moduleColors,
                                savePreservation = TRUE, plot_all_stats = FALSE) {
    multiExpr <- list(ref = list(data = t(exp_ref)), test = list(data = t(exp_test)))
    multiColor <- list(ref = moduleColors)

    # Calculate module preservation
    pres <- WGCNA::modulePreservation(multiExpr, multiColor,
                               referenceNetworks = 1,
                               nPermutations = 200,
                               randomSeed = 1,
                               quickCor = 0,
                               verbose = 3)

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

    # Start the plot
    pdf(file="ModulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
    par(mfrow = c(1,2))
    par(mar = c(4.5,4.5,2.5,1))
    for (p in 1:2)
    {
        min <- min(plotData[, p], na.rm = TRUE)
        max <- max(plotData[, p], na.rm = TRUE)

        # Adjust plotting ranges appropriately
        if (p==2)
        {
            if (min > -max/10) min = -max/10
            ylim <- c(min - 0.1 * (max-min), max + 0.1 * (max-min))
        } else
            ylim <- c(max + 0.1 * (max-min), min - 0.1 * (max-min))
        plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
             main = mains[p],
             cex = 2.4,
             ylab = mains[p], xlab = "Module size", log = "x",
             ylim = ylim,
             xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
        labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);

        # For Zsummary, add threshold lines
        if (p==2) {
            abline(h=0)
            abline(h=2, col = "blue", lty = 2)
            abline(h=10, col = "darkgreen", lty = 2)
        } }

    dev.off()

    if(plot_all_stats == TRUE) {
        # Re-initialize module color labels and sizes
        modColors <- rownames(statsZ)
        moduleSizes <- mp$quality$Z[[ref]][[test]][, 1];

        # Exclude improper modules
        plotMods <- !(modColors %in% c("grey", "gold"));

        # Create numeric labels for each module
        labs <- match(modColors[plotMods], standardColors(50));

        # Start the plot: open a suitably sized graphical window and set sectioning and margins.
        pdf(file = "all_Zstats_combined.pdf", onefile = TRUE)
        par(mfrow = c(3,5))
        par(mar = c(3,3,2,1))
        par(mgp = c(1.6, 0.4, 0));
        # Plot each Z statistic in a separate plot.
        for (s in 1:ncol(statsZ))
        {
            min = min(statsZ[plotMods, s], na.rm = TRUE);
            max = max(statsZ[plotMods, s], na.rm = TRUE);
            if (min > -max/5) min = -max/5
            plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
                 main = colnames(statsZ)[s],
                 cex = 1.7,
                 ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
                 ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
                 xlim = c(20, 1000))
            labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 0.7, offs = 0.04);
            abline(h=0)
            abline(h=2, col = "blue", lty = 2)
            abline(h=10, col = "darkgreen", lty = 2)
        }
        dev.off()
    }
    return(pres)
}

