#' Simple PCA plot for TalusDataSet
#'
#' This simple PC1 vs PC2 PCA helps to check for batch effects and the like.
#'
#' @param se
#' @param top_n
#' @param interest_group
#' @import ggplot2
#' @author Chao-Jen Wong
#' @export

plot_pca <- function(se, top_n = 500, color_by) {

  if (is.list(se)) {
    pcs <- map_dfr(se, function(object) {
      .get_pcs(object, top_n)
    }, .id = "source")

    if (!all(color_by %in% names(pcs))) {
      stop("the argument 'color_by' should specify columns of colData(se)")
    }

    ggplot(pcs, aes_string(x="PC1", y="PC2", color=color_by)) +
      geom_point(size=2, alpha=0.8) +
      facet_wrap(~ source, nrow=2) +
      theme_bw() +
      theme(legend.position    = c(0.7, 0.2))
  }
  else {
    pcs <- .get_pcs(se, top_n)
    ggplot(pcs, aes_string(x="PC1", y="PC2", color=color_by)) +
      geom_point(size=2, alpha=0.8) +
      theme_bw()
  }


}

.get_pcs <- function(object, top_n) {
  # wrap plotPCA from DESeq2
  # calculate the variance for each gene
  rv <- rowVars(assay(object))

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(top_n, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  # assembly the data for the plot
  df <- data.frame(PC1=pca$x[, 'PC1'],
                   PC2=pca$x[, 'PC2'],
                   name=colnames(object), colData(object))

}
