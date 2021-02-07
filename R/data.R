
#' Soybean seed gene expression data from SRA study PRJNA197251
#'
#' Filtered expression data in transcripts per million (TPM) from Danzer et al., 2015.
#' Genes with TPM values < 2 in any of samples were removed to reduce package size.
#' The expression data and associated sample metadata are stored in a SummarizedExperiment
#' object with 10381 genes and 51 samples.
#'
#' @name se.seed
#' @format An object of class \code{SummarizedExperiment}
#' @references Danzer, J., Mellott, E., Bui, A. Q., Le, B. H., Martin, P., Hashimoto, M., ... & Goldberg, R. B. (2015).
#' Down-regulating the expression of 53 soybean transcription factor genes uncovers a role for SPEECHLESS in
#' initiating stomatal cell lineages during embryo development. Plant Physiology, 168(3), 1025-1035.
#' @examples
#' data(se.seed)
"se.seed"


#' Soybean Interpro annotation
#'
#' Interpro protein domain annotation retrieved from the PLAZA 4.0 database.
#' Only genes included in \code{se.seed} are present in this subset.
#'
#' @name soybean_interpro
#' @format A 2-column data frame containing gene IDs and their associated Intepro annotations.
#' @references
#' Van Bel, M., Diels, T., Vancaester, E., Kreft, L., Botzki, A., Van de Peer, Y., ... & Vandepoele, K. (2018). PLAZA 4.0: an integrative resource for functional, evolutionary and comparative plant genomics. Nucleic acids research, 46(D1), D1190-D1196.
#' @examples
#' data(soybean_interpro)
"soybean_interpro"
