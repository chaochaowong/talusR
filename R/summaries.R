#' Summary of Talus limma results
#'
#' Summarize Talus limma results returned by \code{talus_limma}
#' @param object a \code{data.frame} or a list of \code{data.frame} containing limma's \code{lmFit} results.
#' @param alpha  the significance adjust p-vluae cutoff for independent filtering. Default to 0.05.
#' @param lfc_threshold a non-negative values of significance log fold change threshold for independent filtering. Default to NULL.
#'
#' @export
summary_limma <- function(object, alpha = 0.05, lfc_threshold = NULL) {

}
