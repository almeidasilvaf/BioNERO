
#' Soybean seed gene expression data from SRA study PRJNA197251
#'
#' Filtered expression data in transcripts per million (TPM) from Danzer et al., 2015.
#' Genes with TPM values < 1 in any of samples were removed to reduce package size.
#' The expression data and associated sample metadata are stored in a SummarizedExperiment
#' object with 13278 genes and 51 samples.
#'
#' @name se.seed
#' @format An object of class \code{SummarizedExperiment}
#' @references Danzer, J., Mellott, E., Bui, A. Q., Le, B. H., Martin, P., Hashimoto, M., ... & Goldberg, R. B. (2015).
#' Down-regulating the expression of 53 soybean transcription factor genes uncovers a role for SPEECHLESS in
#' initiating stomatal cell lineages during embryo development. Plant Physiology, 168(3), 1025-1035.
#'
#'
"se.seed"
