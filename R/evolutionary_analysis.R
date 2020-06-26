
#' Boolean expression to check if gene or gene set is singleton or not
#'
#' For a given gene or group of genes, this function returns logical values indicating whether they are singletons or not.
#'
#' @param genes Character containing gene or group of genes to be evaluated.
#' @param dup_pairs Data frame containing all duplicate pairs. Column 1 and column 2 must contain gene 1 and gene 2 from the pair, respectively. Addiotional columns are interpreted as pairs info (e.g. mode of duplication).
#'
#' @return Vector of logical values indicating if gene or group of genes is singleton or not.
#'
#' @seealso \code{is_duplicated}
#' @author Fabricio Almeida-Silva
#' @rdname is_singleton
#' @export

is_singleton <- function(genes, dup_pairs) {
    dup_genes <- c(as.character(dup_pairs[,1]), as.character(dup_pairs[,2]))
    eval <- ifelse(genes %in% dup_genes, FALSE, TRUE)
    return(eval)
}

#' Boolean expression to check if gene or gene set is duplicated in the genome or not
#'
#' For a given gene or group of genes, this function returns logical values indicating whether they have copies in the genome or not.
#'
#' @param genes Character containing gene or group of genes to be evaluated.
#' @param dup_pairs Data frame containing all duplicate pairs. Column 1 and column 2 must contain gene 1 and gene 2 from the pair, respectively. Addiotional columns are interpreted as pairs info (e.g. mode of duplication).
#'
#' @return Vector of logical values indicating if gene or group of genes has copies in the genome or not.
#'
#' @seealso \code{is_singleton}
#' @author Fabricio Almeida-Silva
#' @rdname is_duplicated
#' @export

is_duplicated <- function(genes, dup_pairs) {
    dup_genes <- c(as.character(dup_pairs[,1]), as.character(dup_pairs[,2]))
    eval2 <- ifelse(genes %in% dup_genes, TRUE, FALSE)
    return(eval2)
}


#' Get homologous genes of given gene or gene set
#'
#' @param genes Character vector containing gene or set of genes of which homologs will be extracted.
#' @param dup_pairs Data frame containing all ortholog pairs (different species) or paralog pairs (same species). Column 1 and column 2 must contain gene 1 and gene 2 from the pair, respectively. Addiotional columns are interpreted as pairs info (e.g. mode of duplication for paralog pairs).
#'
#' @return List containing homologs of each gene in the input gene set.
#' @author Fabricio Almeida-Silva
#' @export
#' @rdname get_homologs

get_homologs <- function(genes, dup_pairs) {
    filtered_dup <- dup_pairs[dup_pairs[,1] %in% genes | dup_pairs[,2] %in% genes, ]
    genes_and_paralogs <- split(filtered_dup[,1], filtered_dup[,2])
    paralogs <- lapply(genes_and_paralogs, function(x) unique(x[x %in% genes]))
    paralogs <- paralogs[lapply(paralogs, length) > 0]
    return(paralogs)
}
