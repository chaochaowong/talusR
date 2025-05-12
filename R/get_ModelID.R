#' Get ModelID from a stripped cell line name
#'
#' Retrieve the `ModelID` that corresponds to a stripped cell line name in the DepMap model database.
#'
#' @param cell_line_name a characeter string of a cell line name matching to an entry in the "StrippedCellLineName" column in the DepMap's Model database.
#' @return a character string correponding to \code{cell_line_name}.
#' @author Chao-Jen Wong
#' @examples
#' get_ModelID('THP1')
#'
#' @import dplyr
#' @export


get_ModelID <- function(cell_line_name) {
  data(model)
  model_id <- model %>%
    dplyr::filter(str_detect(StrippedCellLineName, cell_line_name)) %>%
    pull(ModelID)

  if (length(model_id) > 0)
    return(model_id)

  if (identical(model_id, character(0)))
    return('ModelID not found.')
}
