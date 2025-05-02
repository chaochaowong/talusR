#'
#'
#' @import limma
#' @import SummarizedExperiment
#' @importFrom tibble rownames_to_column
#' @author Chao-Jen Wong
#' @export
talus_limma <- function(se, design = ~0 + Tx) {
  require(limma)
  # check design formula: does Tx exist, is it a factor
  # display a level and control vs. contrasts

  if (is.list(se)) {
    res <- lapply(se, .wrap_limma, design)
    names(res) <- names(se)

  }
  else
    res <- .wrap_limma(se, design)

  return(res)
}

.wrap_limma <- function(object, design) {
  # assemble row_data to be join to the output of topTable()
  row_data <- as.data.frame(rowData(object)) %>%
    tibble::rownames_to_column(var='id')

  # assemble contrast model matrix
  col_data <- colData(object)
  assay <- assay(object)
  mm <- model.matrix(design, data = col_data)
  fit <- lmFit(assay, mm)
  control <- colnames(mm)[1]
  contrasts <- colnames(mm)[-1]

  # 5. Build a named vector of contrast expressions
  #    e.g. "doxo - DMSO", "alca - DMSO"
  contrast_exprs <- paste0(contrasts, " - ", control)
  names(contrast_exprs) <- paste0(contrasts, "_vs_", control)

  # 6. Create the contrast matrix
  #    you can pass your design matrix directly to levels=
  contrast_matrix <- do.call(makeContrasts,
                         c(as.list(contrast_exprs),
                           list(levels = mm)))

  # 7. Fit contrasts & compute moderated statistics
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  # make data.frame
  top_table <- lapply(names(contrast_exprs), function(conts) {
    topTable(fit2, coef = conts, number = nrow(assay),
             sort.by = 'P') %>%
      tibble::rownames_to_column(var='id') %>%
      dplyr::left_join(row_data, by = 'id')
  })
  names(top_table) <- names(contrast_exprs)

  return(top_table)
}
