#' Simple point plot to show per-protein aboundance
#'
#' Given a protein ID, plot the MS intensity signal
#'
#' @param se SummarizedExperience instance (SE) or list of SE containing the transformed protein abundance.
#' @param protein_id Protein ID
#' @param group_by Group by a factor in the meta file
#'
#' @importFrom purrr map_dfr
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_detect
#' @importFrom dplyr filter
#' @importFrom tidyr gather
#' @import ggplot2
#'
#' @author Chao-Jen Wong
#'
#' @export
#'
per_protein_abun <- function(se, protein_id, category_by) {
  # check if protein_id exist as part of the rownames
  require(ggplot2)
  require(tibble)
  require(tidyr)

  if (is.list(se)) {
    #if( protein_id %in% rownames(se[[1]]))
    # stop(sprintf("Protein ID '%s' not found.", protein_id))
    if (! rlang::has_name(colData(se[[1]]), category_by)) {
      stop(sprintf("Column '%s' not found in %s", which_run, meta))
    }

    df <- map_dfr(se, .extract_per_protein_abun, protein_id,
                  .id = "source")
    ggplot(df, aes_string(x = category_by, y = 'abundance')) +
      geom_point(size=2, alpha = 0.8, color='steelblue') +
      theme_bw() +
      facet_wrap(~source, ncol = 1, scales = 'free_y') +
      labs(title = protein_id) +
      theme(
        panel.grid.minor = element_blank()
      )
  }
  else {
    if (!protein_id %in% rownames(se))
      stop(sprintf("Protein ID '%s' not found.", protein_id))

    df <- .extract_per_protein_abun(se, protein_id)
    ggplot(df, aes_string(x = category_by, y = 'abundance')) +
      geom_point(size=2, alpha = 0.8, color='steelblue') +
      theme_bw() +
      labs(title = protein_id) +
    theme(
      panel.grid.minor = element_blank()
    )
  }
}

.extract_per_protein_abun <- function(object, protein_id) {
  col_data <- as.data.frame(colData(object)) %>%
    tibble::rownames_to_column(var='sample')

  df <- as.data.frame(assay(object)) %>%
    tibble::rownames_to_column(var='id') %>%
    dplyr::filter(str_detect(id, protein_id)) %>%
    tidyr::gather(sample, abundance, -id) %>%
    dplyr::left_join(col_data, by = 'sample')

}
