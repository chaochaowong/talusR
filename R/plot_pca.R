#' Simple PCA plot for TalusDataSet and TalusDataSetList
#'
#' This simple PC1 vs PC2 PCA helps to check for batch effects and the like.
#'
#' @param tds a \code{TalusDataSet} or \code{TalusDataSetList} object storing Talus data
#' @param top_n number of top proteins to use for principal components, selected by highest row variance
#' @param color_by a character string in colData to use for coloring
#' @return
#'   a \code{ggplot} object for \code{TalusDataSet} or \code{TalusDataSetList}
#' @rdname plot_pca
#' @export
setGeneric("plot_pca",
           function(object, ...) standardGeneric("plot_pca"))

#' @rdname plot_pca
#' @aliases plot_pca,TalusDataSetList
#' @exportMethod plot_pca
setMethod("plot_pca", "TalusDataSetList",
          function(object,
                   top_n = 500,
                   color_by) {
    # get pcs
    pcs <- map_dfr(object, function(x) {
      .get_pcs(x, top_n)
    }, .id = "source")

    # make sure color_by column exists
    if (!all(color_by %in% names(pcs))) {
      stop("the argument 'color_by' should specify columns of colData(object)")
    }

    # plot pca
    ggplot(pcs, aes_string(x = "PC1", y = "PC2", color = color_by)) +
      geom_point(size = 2, alpha = 0.8) +
      facet_wrap(~source, nrow = 2, scales = "free") +
      theme_bw() +
      theme(legend.position = c(0.7, 0.2))
  }
)

#' @rdname plot_pca
#' @aliases plot_pca,TalusDataSet
#' @exportMethod plot_pca
setMethod("plot_pca", "TalusDataSet",
          function(object,
                   top_n = 500,
                   color_by) {

    pcs <- .get_pcs(object, top_n)

    if (!all(color_by %in% names(pcs))) {
      stop("the argument 'color_by' should specify columns of colData(object)")
    }

    ggplot(pcs, aes_string(x = "PC1", y = "PC2", color = color_by)) +
      geom_point(size = 2, alpha = 0.8) +
      theme_bw()
  }
)


.get_pcs <- function(object, top_n) {
  # wrap plotPCA from DESeq2
  # calculate the variance for each gene
  rv <- rowVars(assay(object))

  # select the ntop genes by variance
  select <- order(rv, decreasing = TRUE)[seq_len(min(top_n, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select, ]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)

  # assembly the data for the plot
  df <- data.frame(
    PC1 = pca$x[, "PC1"],
    PC2 = pca$x[, "PC2"],
    name = colnames(object), colData(object),
    percent_var = paste(c("PC1:", "PC2:"),
      paste0(round(percentVar * 100)[1:2], "%"),
      collapse = "; "
    )
  )
  return(df)
}
