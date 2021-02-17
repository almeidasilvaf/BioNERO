
#' Maize gene expression data from Shin et al., 2021.
#'
#' Filtered expression data in transcripts per million (TPM) from Shin et al., 2021.
#' Genes with TPM values <5 in more than 60% of the samples were removed to reduce package size.
#' The expression data and associated sample metadata are stored in a SummarizedExperiment object.
#'
#' @name zma.se
#' @format An object of class \code{SummarizedExperiment}
#' @references Shin, J., Marx, H., Richards, A., Vaneechoutte, D., Jayaraman, D., Maeda, J., ... & Roy, S. (2021). A network-based comparative framework to study conservation and divergence of proteomes in plant phylogenies. Nucleic Acids Research, 49(1), e3-e3.
#' @examples
#' data(zma.se)
"zma.se"


#' Rice gene expression data from Shin et al., 2021.
#'
#' Filtered expression data in transcripts per million (TPM) from Shin et al., 2021.
#' Genes with TPM values <5 in more than 60% of the samples were removed to reduce package size.
#' The expression data and associated sample metadata are stored in a SummarizedExperiment object.
#'
#' @name osa.se
#' @format An object of class \code{SummarizedExperiment}
#' @references Shin, J., Marx, H., Richards, A., Vaneechoutte, D., Jayaraman, D., Maeda, J., ... & Roy, S. (2021). A network-based comparative framework to study conservation and divergence of proteomes in plant phylogenies. Nucleic Acids Research, 49(1), e3-e3.
#' @examples
#' data(osa.se)
"osa.se"


#' Filtered maize gene expression data from Shin et al., 2021.
#'
#' Filtered expression data in transcripts per million (TPM) from Shin et al., 2021.
#' This is the same data set described in \code{zma.se}, but it only contains the top 500 genes with the highest variances.
#' This data set was created to be used in unit tests and examples.
#'
#' @name filt.se
#' @format An object of class \code{SummarizedExperiment}
#' @references Shin, J., Marx, H., Richards, A., Vaneechoutte, D., Jayaraman, D., Maeda, J., ... & Roy, S. (2021). A network-based comparative framework to study conservation and divergence of proteomes in plant phylogenies. Nucleic Acids Research, 49(1), e3-e3.
#' @examples
#' data(filt.se)
"filt.se"


#' Maize Interpro annotation
#'
#' Interpro protein domain annotation retrieved from the PLAZA Monocots 4.0 database.
#' Only genes included in \code{zma.se} are present in this subset.
#'
#' @name zma.interpro
#' @format A 2-column data frame containing gene IDs and their associated Interpro annotations.
#' @references
#' Van Bel, M., Diels, T., Vancaester, E., Kreft, L., Botzki, A., Van de Peer, Y., ... & Vandepoele, K. (2018). PLAZA 4.0: an integrative resource for functional, evolutionary and comparative plant genomics. Nucleic acids research, 46(D1), D1190-D1196.
#' @examples
#' data(zma.interpro)
"zma.interpro"


