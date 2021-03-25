
#' Logical expression to check if gene or gene set is singleton or not
#'
#' @param genes Character containing gene or group of genes to be evaluated.
#' @param og Data frame of 3 columns corresponding to orthogroup, species ID,
#' and gene ID, respectively.
#' @return Vector of logical values indicating if gene or group of genes
#' is singleton or not.
#'
#' @seealso \code{is_duplicated}
#' @author Fabricio Almeida-Silva
#' @rdname is_singleton
#' @export
#' @examples
#' data(og.zma.osa)
#' data(filt.se)
#' genes <- tail(rownames(filt.se), n = 100)
#' is_singleton(genes, og.zma.osa)
is_singleton <- function(genes, og) {
    s <- lengths(split(og[og[,3] %in% genes, 3], genes)) == 1
    return(s)
}

