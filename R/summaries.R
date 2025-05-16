#' Summary of Talus differential results
#'
#' Summarize Talus differential results returned by either by \code{talus_limma} or \code{talus_row_t_welch}.
#' @param objecta a list of \code{data.frame} containing limma's \code{lmFit} results of all contrasts.
#' @param alpha the significance adjust p-vluae cutoff for independent filtering. Default to 0.05.
#' @param lfc_threshold a non-negative values of significance log fold change threshold for independent filtering. Default to NULL.
#'
#' @importFrom purrr map
#' @importFrom dplyr filter
#'
#' @export
summary_talus <- function(object,
                          alpha = NULL,
                          lfc_threshold = NULL) {
  if (is.null(alpha) && is.null(lfc_threshold)) {
    stop(
      "At least one of 'alpha' or 'lfc_threshold' must be provided.",
      call. = FALSE
    )
  }

  # if (is.list(object)) {
  # object should be a list
  # get the filtered data.frame
  # display message
  obj_filt <- map2(
    object, names(object), display_summary_talus,
    alpha, lfc_threshold
  )

  # }


  return(obj_filt)
}




display_summary_talus <- function(res,
                                  contrast_name,
                                  alpha,
                                  lfc_threshold) {
  # must check if the input is a data.frame containing either adj.P.Val or logFC columns
  # if (!is.data.frame(res))
  #  stop('')
  filt <- res
  if (!is.null(alpha)) {
    filt <- filt %>%
      dplyr::filter(adj.P.Val < alpha)
  }

  if (!is.null(lfc_threshold)) {
    filt <- filt %>%
      dplyr::filter(abs(logFC) > abs(lfc_threshold))
  }

  total_features <- nrow(res)
  n_up <- sum(filt$logFC > 0, na.rm = TRUE)
  n_down <- sum(filt$logFC < 0, na.rm = TRUE)
  n_zero <- sum(filt$logFC == 0, na.rm = TRUE)

  # print a summary
  cat(contrast_name, ":\n\n")
  cat(sprintf(
    "Threshold: alpha < %s; |lfc_threshod| > %s\n",
    if (is.null(alpha)) "NULL" else alpha,
    if (is.null(lfc_threshold)) "NULL" else lfc_threshold
  ))
  cat(
    sprintf(
      "Total features: %d\nUp-regulated:   %d\n",
      total_features, n_up
    ),
    sprintf("Down-regulated: %d\n", n_down),
    if (n_zero > 0) sprintf("No change:   %d\n", n_zero) else "",
    sep = ""
  )

  cat("\n")
  return(filt)
}
