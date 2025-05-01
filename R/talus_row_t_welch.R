#' Per-protien t-Weltch testing
#'
#' A wrapping of
#' SummarizedExperiment instances with three layers: nuc, cyto, and plasm.
#'
#' @param se A \code{SummarizedExperiment} instance with log-transform assays
#' @param design A formula to create a model matrix. Default to ~0 + Tx, where \code{Tx} is a factor from colData of \code{se}
#'
#' @return Orginal \code{SummarizedExperiment} objects with t-statistics append to rowData of \code{se}
#' @author Chao-Jen Wong
#' @importFrom matrixTests row_t_welch
#' @export
talus_row_t_welch <- function(se, design = ~ 0 + Tx) {
  require(matrixTests)
  # check design formula: does Tx exist, is it a factor
  # display a level and control vs. contrasts

  if (is.list(se)) {
    res <- lapply(se, .wrap_row_t_welch, design)
    names(res) <- names(se)
  }
  else
    res <- .wrap_row_t_welch(se, design)

  return(res)
}

.wrap_row_t_welch <- function(object,
                              design) {
  col_data <- colData(ojbect)
  assay <- assay(object)
  mm <- model.matrix(design, data = col_data)

  control <- colnames(mm)[1]
  contrast <- colnames(mm)[-1]
  control_col <- as.logical(mm[, control, drop=TRUE])

  res <- map_dfc(contrast, function(ct) {
    contrast_col <- as.logical(mm[, ct, drop=TRUE])
    res <- row_t_welch(x = assay[, contrast_col],
                       y = assay[, control_col])
    # tidy up
    res %>%
      dplyr::select(mean.x, mean.diff, statistic, pvalue) %>%
      dplyr::mutate(padj = p.adjust(pvalue, method = 'BH')) %>%
      dplyr::rename(avg_log2 = mean.x,
                    lfc = mean.diff) %>%
      rename_with(~ paste0(ct, '_', .))
  }) %>%
    dplyr::mutate(control_mean  = rowMeans(assay[, control_col]),
                  .before = 1) %>%
    dplyr::rename(!!paste0(control, '_avg_log2') := control_mean)

  return(res)
}
