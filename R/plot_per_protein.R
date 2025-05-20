#' Simple point plot to show per-protein aboundance
#'
#' Given a protein ID, plot the MS intensity signal
#'
#' @param object a \code{TalusDataSet} or \code{TalusDataSetList} instance containing the transformed protein abundance.
#' @param protein_id Protein ID
#' @param group_by group by a column (factor) in the meta file when plotting
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
#' @rdname plot_per_protein
setGeneric("plot_per_protein",
           function(object, protein_id, group_by)
             standardGeneric("plot_per_protein"))

#' @rdname plot_per_protein
#' @aliases plot_per_protein,TalusDataSetList
#' @exportMethod plot_per_protein
setMethod("plot_per_protein", "TalusDataSetList",
   function(object,
            protein_id,
            group_by) {

    if (!protein_id %in% rownames(object[[1]]))
      stop(sprintf("Protein ID '%s' not found.", protein_id))

    col_data <- colData(object[[1]])
    if (!rlang::has_name(colData(object[[1]]), group_by)) {
      stop(sprintf("Column '%s' not found in the colData of the object", group_by))
    }

    df <- map_dfr(object, .extract_per_protein, protein_id,
      .id = "source"
    )
    ggplot(df, aes_string(x = group_by, y = "abundance")) +
      geom_point(size = 2, alpha = 0.8, color = "steelblue") +
      theme_bw() +
      facet_wrap(~source, ncol = 1, scales = "free_y") +
      labs(title = protein_id) +
      theme(
        panel.grid.minor = element_blank()
      )
   }
)

#' @rdname plot_per_protein
#' @aliases plot_per_protein,TalusDataSet
#' @exportMethod plot_per_protein
setMethod("plot_per_protein", "TalusDataSet",
  function(object,
           protein_id,
           group_by) {

    if (protein_id %in% rownames(object))
      stop(sprintf("Protein ID '%s' not found.", protein_id))

    col_data <- colData(object)
    if (!rlang::has_name(col_data, group_by)) {
      stop(sprintf("Column '%s' not found in the colData of the object", group_by))
    }


    df <- .extract_per_protein(object, protein_id)
    ggplot(df, aes_string(x = group_by, y = "abundance")) +
      geom_point(size = 2, alpha = 0.8, color = "steelblue") +
      theme_bw() +
      labs(title = protein_id) +
      theme(
        panel.grid.minor = element_blank()
      )
  }
)

.extract_per_protein <- function(object, protein_id) {
  col_data <- as.data.frame(colData(object)) %>%
    tibble::rownames_to_column(var = "sample")

  df <- as.data.frame(assay(object)) %>%
    tibble::rownames_to_column(var = "id") %>%
    dplyr::filter(str_detect(id, protein_id)) %>%
    tidyr::gather(sample, abundance, -id) %>%
    dplyr::left_join(col_data, by = "sample")
}
