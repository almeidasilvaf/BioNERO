#' Identify consensus modules across multiple data sets
#'
#' This function identifies modules that are conserved across multiple data sets. It can be used for evolutionary analyses and for testing robustness of coexpression modules in data sets from different sources, for instance.
#'
#' @param exp_list A list containing the expression data frames. For each expression data frame, row names correspond to gene IDs and column names correspond to sample names. The list can be created by using \code{list(exp1, exp2, ..., expn)}.
#' @param setLabels Character vector containing labels for each expression set.
#' @param metadata Data frame containing sample information. The first column must contain sample names and the second column must contain sample information.
#' @param cor_method Correlation method used for network reconstruction. One of "spearman" (default), "biweight", or "pearson".
#' @param net_type Network type. One of "signed hybrid" (default), "signed" or "unsigned".
#' @param rsquared Minimum R squared to consider the network similar to a scale-free topology. Default is 0.8.
#' @param module_merging_threshold Cut-off for merging similar modules based on the dissimilarity of their eigengenes. Default is 0.2, indicating that modules that have correlation >= 0.8 will be merged into one.
#'
#' @return A list containing 4 elements: \describe{
#'   \item{consModules}{Consensus module assignments}
#'   \item{consMEs}{Consensus module eigengenes}
#'   \item{exprSize}{Description of the multi-set object returned by the function \code{WGCNA::checkSets}}
#'   \item{sampleInfo}{Metadata for each expression set}
#' }
#'
#' @seealso
#'  \code{\link[WGCNA]{checkSets}},\code{\link[WGCNA]{pickSoftThreshold}},\code{\link[WGCNA]{adjacency}},\code{\link[WGCNA]{TOMsimilarity}},\code{\link[WGCNA]{pquantile}},\code{\link[WGCNA]{labels2colors}},\code{\link[WGCNA]{multiSetMEs}},\code{\link[WGCNA]{consensusMEDissimilarity}},\code{\link[WGCNA]{mergeCloseModules}},\code{\link[WGCNA]{plotDendroAndColors}}
#'  \code{\link[pals]{discrete}}
#'  \code{\link[dynamicTreeCut]{cutreeDynamic}}
#' @rdname consensus_modules
#' @export
#' @importFrom WGCNA checkSets collectGarbage pickSoftThreshold addGrid adjacency TOMsimilarity pquantile labels2colors multiSetMEs consensusMEDissimilarity mergeCloseModules plotDendroAndColors
#' @importFrom pals stepped
#' @importFrom dynamicTreeCut cutreeDynamic
consensus_modules <- function(exp_list, setLabels = NULL, metadata,
                              cor_method = "spearman",
                              net_type = "signed hybrid", rsquared = 0.8,
                              module_merging_threshold = 0.2) {
    nSets <- length(exp_list)

    # Build multi-set object
    multiExp <- vector(mode = "list", length = nSets)
    for (set in 1:nSets) {
        multiExp[[set]] <- list(data = as.data.frame(t(exp_list[[set]])))
    }

    # Check if the data has the correct format for downstream analysis
    expSize <- WGCNA::checkSets(multiExp)

    # Build multi-set object with metadata
    sampleinfo <- vector(mode = "list", length = nSets)
    for (set in 1:nSets) {
        setSamples <- rownames(multiExp[[set]]$data)
        traitRows <- match(setSamples, metadata[,1])
        sampleinfo[[set]] <- list(data = metadata[traitRows, -1])
        rownames(sampleinfo[[set]]$data) <- metadata[traitRows, 1]
    }

    WGCNA::collectGarbage()

    # Define dimensions
    nGenes <- expSize$nGenes
    nSamples <- expSize$nSamples

    # Choose a set of soft-thresholding powers
    powers <- c(seq(4, 10, by = 1), seq(12, 20, by = 2))
    sft <- vector(mode = "list", length = nSets)
    sft_power <- vector(mode = "list", length = nSets)
    powerTables <- vector(mode = "list", length = nSets)
    for(set in 1:nSets) {
        if(cor_method == "pearson") {
            sft[[set]] = WGCNA::pickSoftThreshold(multiExp[[set]]$data,
                                                  networkType = net_type, powerVector=powers,
                                                  verbose = 2, RsquaredCut = rsquared)
            powerTables[[set]] <- list(data = sft[[set]]$fitIndices)
            sft_power[[set]] <- sft[[set]]$powerEstimate

        } else if(cor_method == "biweight") {
            sft[[set]] <- WGCNA::pickSoftThreshold(multiExp[[set]]$data,
                                                   networkType = net_type, powerVector=powers, RsquaredCut = rsquared,
                                                   corFnc = bicor, corOptions = list(use = 'p', maxPOutliers = 0.05))
            powerTables[[set]] <- list(data = sft[[set]]$fitIndices)
            sft_power[[set]] <- sft[[set]]$powerEstimate

        } else if (cor_method == "spearman") {
            sft[[set]] <- WGCNA::pickSoftThreshold(multiExp[[set]]$data,
                                                   networkType = net_type, powerVector=powers, RsquaredCut = rsquared,
                                                   corOptions = list(use = 'p', method = "spearman"))
            powerTables[[set]] <- list(data = sft[[set]]$fitIndices)
            sft_power[[set]] <- sft[[set]]$powerEstimate

        } else {
            print("Please, specify a correlation method (one of 'spearman', 'pearson' or 'biweight').")
        }
    }
    WGCNA::collectGarbage()

    # Plot the results:
    colors = pals::stepped(n = nSets)
    plotCols = c(2,5,6,7)
    colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
                 "Max connectivity")

    # Get the minima and maxima of the plotted points
    ylim = matrix(NA, nrow = 2, ncol = 4)
    for (set in 1:nSets)
    {
        for (col in 1:length(plotCols))
        {
            ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
            ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
        }
    }
    # Plot the quantities in the chosen columns vs. the soft thresholding power
    pdf(file = "SFT_fit_consensus.pdf", width = 8, height = 6)
    par(mfcol = c(2,2))
    par(mar = c(4.2, 4.2 , 2.2, 0.5))
    cex1 = 0.7
    for (col in 1:length(plotCols)) for (set in 1:nSets)
    {
        if (set==1)
        {
            plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
                 xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
                 main = colNames[col])
            WGCNA::addGrid()
        }
        if (col==1)
        {
            text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
                 labels=powers,cex=cex1,col=colors[set])
        } else
            text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
                 labels=powers,cex=cex1,col=colors[set])
        if (col==1)
        {
            legend("bottomright", legend = setLabels, col = colors, pch = 20)
        } else
            legend("topright", legend = setLabels, col = colors, pch = 20)
    }
    dev.off()

    # Calculate adjacencies for each individual set
    softPower <- sft_power
    adjacencies <- array(0, dim = c(nSets, nGenes, nGenes))
    print("Calculating adjacency matrix...")
    for (set in 1:nSets) {
        if(cor_method == "pearson") {
            adjacencies[set, , ] <- WGCNA::adjacency(multiExp[[set]]$data, power=softPower[[set]], type=net_type)
        } else if(cor_method == "spearman") {
            adjacencies[set, , ] <- WGCNA::adjacency(multiExp[[set]]$data, power=softPower[[set]], type=net_type,
                                                     corOptions = list(use = "p", method = "spearman"))
        } else if (cor_method == "biweight") {
            adjacencies[set, , ] <- WGCNA::adjacency(multiExp[[set]]$data, power=softPower[[set]], type=net_type,
                                                     corFnc = bicor)
        } else {
            print("Please, specify a correlation method. One of 'spearman', 'pearson' or 'biweight'.")
        }
    }

    # Calculate TOMs in each individual data set
    TOM <- array(0, dim = c(nSets, nGenes, nGenes))
    for (set in 1:nSets) {
        #Calculate TOM from adjacency matrix
        print("Calculating topological overlap matrix (TOM)...")
        if(net_type == "signed hybrid") {
            TOM[set, , ] <- WGCNA::TOMsimilarity(adjacencies[set, , ], TOMType = "signed")
        } else if(net_type == "signed") {
            TOM[set, , ] <- WGCNA::TOMsimilarity(adjacencies[set, , ], TOMType = "signed Nowick")
        } else {
            TOM[set, , ] <- WGCNA::TOMsimilarity(adjacencies[set, , ], TOMType = "unsigned")
        }
    }

    # Scaling TOM to make them comparable
    scaleP <- 0.95
    set.seed(123)

    # Sample sufficiently large number of TOM entries
    nSamples <- as.integer(1 / (1-scaleP) * 1000)

    # Choose the sampled TOM entries
    scaleSample <- sample(nGenes * (nGenes-1) / 2, size = nSamples)
    TOMScalingSamples <- list()

    # TOM values at reference percentile
    scaleQuant <- rep(1, nSets)

    # Scaling powers to equalize TOM values
    scalePowers <- rep(1, nSets)

    # Loop over sets
    for (set in 1:nSets) {

        #Select the sampled TOM entries
        TOMScalingSamples[[set]] <- as.dist(TOM[set, , ])[scaleSample]

        #Calculate the 95th percentile
        scaleQuant[set] <- quantile(TOMScalingSamples[[set]],
                                    probs = scaleP, type = 8)

        # Scale one TOM
        if (set > 1) {
            scalePowers[set] <- log(scaleQuant[1]) / log(scaleQuant[set])
            TOM[set, , ] <- TOM[set, , ]^scalePowers[set]
        }
    }

    # Plot the unscaled TOM VS. the scaled TOM
    scaledTOMSamples <- list()
    for (set in 1:nSets) {
        scaledTOMSamples[[set]] <- TOMScalingSamples[[set]]^scalePowers[set]
    }


    if (nSets == 2) {
        pdf(file = "TOMScaling-QQplot.pdf", width=6, height=6)
        qqUnscaled <- qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]],
                             plot.it=TRUE, cex=0.6,
                             xlab=paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[2]),
                             main = "QQ plot of TOM", pch=20)

        qqScaled <- qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]],
                           plot.it = FALSE)
        points(qqScaled$x, qqScaled$y, col = colors, cex = 0.6, pch = 20)
        abline(a=0, b=1, col = "blue")
        legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch=20,
               col = c("black", "red"))
        dev.off()
    } else {
        print("QQ-plot will not be plotted. More than 2 data sets.")
    }

    # Calculation of consensus TOM
    if (nSets <= 3) {
        consensusTOM <- do.call(pmin, lapply(seq(dim(TOM)[1]), function(i) TOM[i,,]))
    } else {
        consensusTOM <- do.call(function(x) WGCNA::pquantile(x, prob=0.25),
                                lapply(seq(dim(TOM)[1]), function(i) TOM[i,,]))
    }

    # Clustering and module identification
    consTree <- hclust(as.dist(1-consensusTOM), method = "average")
    unmergedLabels <- dynamicTreeCut::cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
                                                    deepSplit = 2, cutHeight = 0.995,
                                                    minClusterSize = 30,
                                                    pamRespectsDendro = FALSE)

    unmergedColors <- WGCNA::labels2colors(unmergedLabels)

    # Merge similar modules
    unmergedMEs <- WGCNA::multiSetMEs(multiExp, colors = NULL, universalColors = unmergedColors)
    consMEDiss <- WGCNA::consensusMEDissimilarity(unmergedMEs)
    consMETree <- hclust(as.dist(consMEDiss), method = "average")

    merge <- WGCNA::mergeCloseModules(multiExp, unmergedLabels,
                                      cutHeight = module_merging_threshold, verbose = 3)

    moduleLabels <- merge$colors
    moduleColors <- WGCNA::labels2colors(moduleLabels)
    consMEs <- merge$newMEs

    # Plot dendrogram with merged colors
    pdf(file = "Consensus_modules_before_and_after_merging.pdf", width=9, height=6)
    WGCNA::plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors),
                               c("Unmerged", "Merged"),
                               dendroLabels = FALSE, hang = 0.03,
                               addGuide = TRUE, guideHang = 0.05)
    dev.off()

    result_list <- list(consModules = moduleColors,
                        consMEs = consMEs,
                        exprSize = expSize,
                        sampleInfo = sampleinfo)

    return(result_list)
}


#' Relate consensus modules to set-specific modules
#'
#' @param setMEs Module eigengenes for the set-specific network. It is the 2nd element of the result list from \code{exp2net}.
#' @param setColors Module assignment for the set-specific network. It is the 5th element of the result list from \code{exp2net}.
#' @param consMEs Module eigengenes for consensus network. It is the 2nd element of the result list from \code{consensus_modules}.
#' @param consColors Module assignment for the consensus network. It is the 1st element of the result list from \code{consensus_modules}.
#' @return Heatmap of relationships between consensus modules and set-specific modules in a PDF file.
#' @seealso
#'  \code{\link[WGCNA]{orderMEs}},\code{\link[WGCNA]{labeledHeatmap}},\code{\link[WGCNA]{blueWhiteRed}}
#' @rdname consensus_set_relationship
#' @export
#' @importFrom WGCNA orderMEs labeledHeatmap blueWhiteRed
consensus_set_relationship <- function(setMEs, setColors, consMEs, consColors) {
    setMEs <- WGCNA::orderMEs(MEs)

    setModuleLabels <- substring(names(setMEs), 3)
    consModuleLabels <- substring(names(consMEs[[1]]$data), 3)

    # Convert the numeric module labels to color labels
    setModules <- labels2colors(as.numeric(setModuleLabels))
    consModules <- labels2colors(as.numeric(consModuleLabels))

    # Numbers of female and consensus modules
    nsetMods <- length(setModules)
    nConsMods <- length(consModules)

    # Initialize tables of p-values and of the corresponding counts
    pTable <- matrix(0, nrow = nsetMods, ncol = nConsMods)
    CountTbl <- matrix(0, nrow = nsetMods, ncol = nConsMods)

    # Execute all pairwaise comparisons
    for (smod in 1:nsetMods)
        for (cmod in 1:nConsMods)
        {
            setMembers <- (setColors == setModules[smod])
            consMembers <- (consColors == consModules[cmod])
            pTable[smod, cmod] <- -log10(fisher.test(femMembers, consMembers, alternative = "greater")$p.value)
            CountTbl[smod, cmod] <- sum(setColors == setModules[smod] & consColors ==
                                            consModules[cmod])
        }

    # Truncate p values smaller than 10^{-50} to 10^{-50}
    pTable[is.infinite(pTable)] <- 1.3 * max(pTable[is.finite(pTable)])
    pTable[pTable>50 ] <- 50

    # Marginal counts (really module sizes)
    setModTotals <- apply(CountTbl, 1, sum)
    consModTotals <- apply(CountTbl, 2, sum)

    # Actual plotting
    pdf(file = "ConsensusVsSet-specificModules.pdf", width = 10, height = 7)
    par(mfrow = c(1,1))
    par(cex = 1.0)
    par(mar = c(8, 10.4, 2.7, 1) + 0.3)

    # Use function labeledHeatmap to produce the color-coded table with all the trimmings
    WGCNA::labeledHeatmap(Matrix = pTable,
                          xLabels = paste(" ", consModules),
                          yLabels = paste(" ", setModules),
                          colorLabels = TRUE,
                          xSymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
                          ySymbols = paste("Set ", setModules, ": ", setModTotals, sep=""),
                          textMatrix = CountTbl,
                          colors = WGCNA::blueWhiteRed(100)[50:100],
                          main = "Correspondence of set-specific and consensus modules",
                          cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE)
    dev.off()

}


#' Correlate set-specific modules and consensus modules to sample information
#'
#' @param consMEs Data frame of consensus module eigengenes returned by \code{consensus_modules}.
#' @param exprSize Object containing information on the multi-set used for inference of consensus modules. Third object of the result list from \code{consensus_modules}.
#' @param sampleInfo Data frame of two columns containing sample names in the first column and sample information in the second column.
#' @param cor_method Correlation method to be used. One of 'spearman' or 'pearson'. Default is 'spearman'.
#' @return Heatmaps of relationships between consensus/set-specific modules and sample information (metadata) in a PDF file.
#'
#' @seealso
#'  \code{\link[WGCNA]{corPvalueFisher}},\code{\link[WGCNA]{labels2colors}},\code{\link[WGCNA]{labeledHeatmap}},\code{\link[WGCNA]{blueWhiteRed}}
#' @rdname consensusmodules_sample_cor
#' @export
#' @importFrom WGCNA corPvalueFisher labels2colors labeledHeatmap blueWhiteRed
consensusmodules_sample_cor <- function(consMEs, exprSize, sampleInfo, cor_method = "spearman") {
    moduleTraitCor <- list()
    moduleTraitPvalue <- list()
    nSets <- exprSize$nSets

    # Calculate the correlations
    for (set in 1:nSets) {
        if(cor_method == "spearman") {
            moduleTraitCor[[set]] = cor(consMEs[[set]]$data, sampleInfo[[set]]$data, use = "p", method = "spearman")
            moduleTraitPvalue[[set]] = WGCNA::corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set])
        } else if(cor_method == "pearson") {
            moduleTraitCor[[set]] = cor(consMEs[[set]]$data, Traits[[set]]$data, use = "p", method = "pearson")
            moduleTraitPvalue[[set]] = WGCNA::corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set])
        }
    }

    MEColors <- WGCNA::labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)))
    MEColorNames <- paste("ME", MEColors, sep="")

    # Plot the module-trait relationship table for sets
    pdf(file = "ModuleSampleRelationships.pdf", width = 10, height = 7, onefile = TRUE)
    for (set in 1:nSets) {
        textMatrix <- paste(signif(moduleTraitCor[[set]], 2), "\n(",
                            signif(moduleTraitPvalue[[set]], 1), ")", sep = "")
        dim(textMatrix) <- dim(moduleTraitCor[[set]])
        par(mar = c(6, 8.8, 3, 2.2))
        WGCNA::labeledHeatmap(Matrix = moduleTraitCor[[set]],
                              xLabels = names(sampleInfo[[set]]$data),
                              yLabels = MEColorNames,
                              ySymbols = MEColorNames,
                              colorLabels = FALSE,
                              colors = blueWhiteRed(50),
                              textMatrix = textMatrix,
                              setStdMargins = FALSE,
                              cex.text = 0.5,
                              zlim = c(-1,1),
                              main = paste("Module-sample relationships"))
    }
    dev.off()

    # Initialize matrices to hold the consensus correlation and p-value
    consensusCor <- matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
    consensusPvalue <- matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))

    if(nSets == 2) {
        negative <- moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0
        consensusCor[negative] <- pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative])
        consensusPvalue[negative] <- pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative])

        positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0
        consensusCor[positive] <- pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive])
        consensusPvalue[positive] <- pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive])
    } else if(nSets == 3) {
        negative <- moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0 & moduleTraitCor[[3]] < 0
        consensusCor[negative] <- pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative],
                                       moduleTraitCor[[3]][negative])
        consensusPvalue[negative] <- pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative],
                                          moduleTraitCor[[3]][negative])

        positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0 & moduleTraitCor[[3]] > 0
        consensusCor[positive] <- pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive],
                                       moduleTraitCor[[3]][positive])
        consensusPvalue[positive] <- pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive],
                                          moduleTraitCor[[3]][positive])
    } else if(nSets == 4) {
        negative <- moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0 & moduleTraitCor[[3]] < 0 &
            moduleTraitCor[[4]] < 0
        consensusCor[negative] <- pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative],
                                       moduleTraitCor[[3]][negative], moduleTraitCor[[4]][negative])
        consensusPvalue[negative] <- pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative],
                                          moduleTraitCor[[3]][negative], moduleTraitCor[[4]][negative])

        positive <- moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0 & moduleTraitCor[[3]] > 0 &
            moduleTraitCor[[4]] > 0
        consensusCor[positive] <- pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive],
                                       moduleTraitCor[[3]][positive], moduleTraitCor[[4]][positive])
        consensusPvalue[positive] <- pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive],
                                          moduleTraitCor[[3]][positive], moduleTraitCor[[4]][positive])
    } else if(nSets == 5) {
        negative <- moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0 & moduleTraitCor[[3]] < 0 &
            moduleTraitCor[[4]] < 0 & moduleTraitCor[[5]] < 0
        consensusCor[negative] <- pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative],
                                       moduleTraitCor[[3]][negative], moduleTraitCor[[4]][negative],
                                       moduleTraitCor[[5]][negative])
        consensusPvalue[negative] <- pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative],
                                          moduleTraitCor[[3]][negative], moduleTraitCor[[4]][negative],
                                          moduleTraitCor[[5]][negative])

        positive <- moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0 & moduleTraitCor[[3]] > 0 &
            moduleTraitCor[[3]] > 0 & moduleTraitCor[[5]] > 0
        consensusCor[positive] <- pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive],
                                       moduleTraitCor[[3]][positive], moduleTraitCor[[4]][positive],
                                       moduleTraitCor[[5]][positive])
        consensusPvalue[positive] <- pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive],
                                          moduleTraitCor[[3]][positive], moduleTraitCor[[4]][positive],
                                          moduleTraitCor[[5]][positive])
    } else if(nSets == 6) {
        negative <- moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0 & moduleTraitCor[[3]] < 0 &
            moduleTraitCor[[4]] < 0 & moduleTraitCor[[5]] < 0 & moduleTraitCor[[6]] < 0
        consensusCor[negative] <- pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative],
                                       moduleTraitCor[[3]][negative], moduleTraitCor[[4]][negative],
                                       moduleTraitCor[[5]][negative], moduleTraitCor[[6]][negative])
        consensusPvalue[negative] <- pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative],
                                          moduleTraitCor[[3]][negative], moduleTraitCor[[4]][negative],
                                          moduleTraitCor[[5]][negative], moduleTraitCor[[6]][negative])

        positive <- moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0 & moduleTraitCor[[3]] > 0 &
            moduleTraitCor[[3]] > 0 & moduleTraitCor[[5]] > 0 & moduleTraitCor[[6]] > 0
        consensusCor[positive] <- pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive],
                                       moduleTraitCor[[3]][positive], moduleTraitCor[[4]][positive],
                                       moduleTraitCor[[5]][positive], moduleTraitCor[[6]][positive])
        consensusPvalue[positive] <- pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive],
                                          moduleTraitCor[[3]][positive], moduleTraitCor[[4]][positive],
                                          moduleTraitCor[[5]][positive], moduleTraitCor[[6]][positive])
    } else if(nSets == 7) {
        negative <- moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0 & moduleTraitCor[[3]] < 0 &
            moduleTraitCor[[4]] < 0 & moduleTraitCor[[5]] < 0 & moduleTraitCor[[6]] < 0 &
            moduleTraitCor[[7]] < 0
        consensusCor[negative] <- pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative],
                                       moduleTraitCor[[3]][negative], moduleTraitCor[[4]][negative],
                                       moduleTraitCor[[5]][negative], moduleTraitCor[[6]][negative],
                                       moduleTraitCor[[7]][negative])
        consensusPvalue[negative] <- pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative],
                                          moduleTraitCor[[3]][negative], moduleTraitCor[[4]][negative],
                                          moduleTraitCor[[5]][negative], moduleTraitCor[[6]][negative],
                                          moduleTraitCor[[7]][negative])

        positive <- moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0 & moduleTraitCor[[3]] > 0 &
            moduleTraitCor[[3]] > 0 & moduleTraitCor[[5]] > 0 & moduleTraitCor[[6]] > 0 &
            moduleTraitCor[[7]] > 0
        consensusCor[positive] <- pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive],
                                       moduleTraitCor[[3]][positive], moduleTraitCor[[4]][positive],
                                       moduleTraitCor[[5]][positive], moduleTraitCor[[6]][positive],
                                       moduleTraitCor[[7]][positive])
        consensusPvalue[positive] <- pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive],
                                          moduleTraitCor[[3]][positive], moduleTraitCor[[4]][positive],
                                          moduleTraitCor[[5]][positive], moduleTraitCor[[6]][positive],
                                          moduleTraitCor[[7]][positive])
    } else if(nSets == 8) {
        negative <- moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0 & moduleTraitCor[[3]] < 0 &
            moduleTraitCor[[4]] < 0 & moduleTraitCor[[5]] < 0 & moduleTraitCor[[6]] < 0 &
            moduleTraitCor[[7]] < 0 & moduleTraitCor[[8]] < 0
        consensusCor[negative] <- pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative],
                                       moduleTraitCor[[3]][negative], moduleTraitCor[[4]][negative],
                                       moduleTraitCor[[5]][negative], moduleTraitCor[[6]][negative],
                                       moduleTraitCor[[7]][negative], moduleTraitCor[[8]][negative])
        consensusPvalue[negative] <- pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative],
                                          moduleTraitCor[[3]][negative], moduleTraitCor[[4]][negative],
                                          moduleTraitCor[[5]][negative], moduleTraitCor[[6]][negative],
                                          moduleTraitCor[[7]][negative], moduleTraitCor[[8]][negative])

        positive <- moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0 & moduleTraitCor[[3]] > 0 &
            moduleTraitCor[[3]] > 0 & moduleTraitCor[[5]] > 0 & moduleTraitCor[[6]] > 0 &
            moduleTraitCor[[7]] > 0 & moduleTraitCor[[8]] > 0
        consensusCor[positive] <- pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive],
                                       moduleTraitCor[[3]][positive], moduleTraitCor[[4]][positive],
                                       moduleTraitCor[[5]][positive], moduleTraitCor[[6]][positive],
                                       moduleTraitCor[[7]][positive], moduleTraitCor[[8]][positive])
        consensusPvalue[positive] <- pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive],
                                          moduleTraitCor[[3]][positive], moduleTraitCor[[4]][positive],
                                          moduleTraitCor[[5]][positive], moduleTraitCor[[6]][positive],
                                          moduleTraitCor[[7]][positive], moduleTraitCor[[8]][positive])
    } else if(nSets == 9) {
        negative <- moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0 & moduleTraitCor[[3]] < 0 &
            moduleTraitCor[[4]] < 0 & moduleTraitCor[[5]] < 0 & moduleTraitCor[[6]] < 0 &
            moduleTraitCor[[7]] < 0 & moduleTraitCor[[8]] < 0 & moduleTraitCor[[9]] < 0
        consensusCor[negative] <- pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative],
                                       moduleTraitCor[[3]][negative], moduleTraitCor[[4]][negative],
                                       moduleTraitCor[[5]][negative], moduleTraitCor[[6]][negative],
                                       moduleTraitCor[[7]][negative], moduleTraitCor[[8]][negative],
                                       moduleTraitCor[[9]][negative])
        consensusPvalue[negative] <- pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative],
                                          moduleTraitCor[[3]][negative], moduleTraitCor[[4]][negative],
                                          moduleTraitCor[[5]][negative], moduleTraitCor[[6]][negative],
                                          moduleTraitCor[[7]][negative], moduleTraitCor[[8]][negative],
                                          moduleTraitCor[[9]][negative])

        positive <- moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0 & moduleTraitCor[[3]] > 0 &
            moduleTraitCor[[3]] > 0 & moduleTraitCor[[5]] > 0 & moduleTraitCor[[6]] > 0 &
            moduleTraitCor[[7]] > 0 & moduleTraitCor[[8]] > 0 & moduleTraitCor[[9]] > 0
        consensusCor[positive] <- pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive],
                                       moduleTraitCor[[3]][positive], moduleTraitCor[[4]][positive],
                                       moduleTraitCor[[5]][positive], moduleTraitCor[[6]][positive],
                                       moduleTraitCor[[7]][positive], moduleTraitCor[[8]][positive],
                                       moduleTraitCor[[9]][positive])
        consensusPvalue[positive] <- pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive],
                                          moduleTraitCor[[3]][positive], moduleTraitCor[[4]][positive],
                                          moduleTraitCor[[5]][positive], moduleTraitCor[[6]][positive],
                                          moduleTraitCor[[7]][positive], moduleTraitCor[[8]][positive],
                                          moduleTraitCor[[9]][positive])
    } else if(nSets == 10) {
        negative <- moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0 & moduleTraitCor[[3]] < 0 &
            moduleTraitCor[[4]] < 0 & moduleTraitCor[[5]] < 0 & moduleTraitCor[[6]] < 0 &
            moduleTraitCor[[7]] < 0 & moduleTraitCor[[8]] < 0 & moduleTraitCor[[9]] < 0 &
            moduleTraitCor[[10]] < 0
        consensusCor[negative] <- pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative],
                                       moduleTraitCor[[3]][negative], moduleTraitCor[[4]][negative],
                                       moduleTraitCor[[5]][negative], moduleTraitCor[[6]][negative],
                                       moduleTraitCor[[7]][negative], moduleTraitCor[[8]][negative],
                                       moduleTraitCor[[9]][negative], moduleTraitCor[[10]][negative])
        consensusPvalue[negative] <- pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative],
                                          moduleTraitCor[[3]][negative], moduleTraitCor[[4]][negative],
                                          moduleTraitCor[[5]][negative], moduleTraitCor[[6]][negative],
                                          moduleTraitCor[[7]][negative], moduleTraitCor[[8]][negative],
                                          moduleTraitCor[[9]][negative], moduleTraitCor[[10]][negative])

        positive <- moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0 & moduleTraitCor[[3]] > 0 &
            moduleTraitCor[[3]] > 0 & moduleTraitCor[[5]] > 0 & moduleTraitCor[[6]] > 0 &
            moduleTraitCor[[7]] > 0 & moduleTraitCor[[8]] > 0 & moduleTraitCor[[9]] > 0 &
            moduleTraitCor[[10]] > 0
        consensusCor[positive] <- pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive],
                                       moduleTraitCor[[3]][positive], moduleTraitCor[[4]][positive],
                                       moduleTraitCor[[5]][positive], moduleTraitCor[[6]][positive],
                                       moduleTraitCor[[7]][positive], moduleTraitCor[[8]][positive],
                                       moduleTraitCor[[9]][positive], moduleTraitCor[[10]][positive])
        consensusPvalue[positive] <- pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive],
                                          moduleTraitCor[[3]][positive], moduleTraitCor[[4]][positive],
                                          moduleTraitCor[[5]][positive], moduleTraitCor[[6]][positive],
                                          moduleTraitCor[[7]][positive], moduleTraitCor[[8]][positive],
                                          moduleTraitCor[[9]][positive], moduleTraitCor[[10]][positive])
    } else {
        stop("Maximum number of expression data sets for consensus module analysis is 10.")
    }

    textMatrix <- paste(signif(consensusCor, 2), "\n(",
                        signif(consensusPvalue, 1), ")", sep = "")
    dim(textMatrix) <- dim(moduleTraitCor[[set]])
    pdf(file = "ModuleTraitRelationships-consensus.pdf", width = 10, height = 7)
    par(mar = c(6, 8.8, 3, 2.2))
    WGCNA::labeledHeatmap(Matrix = consensusCor,
                          xLabels = names(Traits[[set]]$data),
                          yLabels = MEColorNames,
                          ySymbols = MEColorNames,
                          colorLabels = FALSE,
                          colors = WGCNA::blueWhiteRed(50),
                          textMatrix = textMatrix,
                          setStdMargins = FALSE,
                          cex.text = 0.5,
                          zlim = c(-1,1),
                          main = paste("Consensus module--trait relationships across\n",
                                       paste(setLabels, collapse = " and ")))

}


