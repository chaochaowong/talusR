#’ Differential analysis via limma for TalusDataSet objects
#’
#’ Fit linear models and contrasts using \code{limma::lmFit} and \code{limma::eBayes}
#’ on a \code{TalusDataSet} or each element of a \code{TalusDataSetList}.
#' @param object a \code{TalusDataSet} or \code{TalusDataSetList} instance containing protein log-transform assays
#' @param desgin a formula to create a model matrix. Default to ~0 + Tx, where \code{Tx} is a factor from colData of \code{se}
#' @return a \code{TalusResult} if \code{object} is a \code{TalusDataSet}; or a \code{TalusResultList} if \code{TalusDataSetList}.
#'
#' @import limma
#' @importFrom tibble rownames_to_column
#' @author Chao-Jen Wong
#' @name talus_limma
#' @rdname talus_limma
#' @export
setGeneric("talus_limma",
           function(object, design = ~ 0 + Tx)
             standardGeneric("talus_limma")
)

#' @rdname talus_limma
#' @aliases talus_limma,TalusDataSetList
#' @exportMethod talus_limma
setMethod("talus_limma", signature(object = "TalusDataSetList"),
  function(object, design = ~ 0 + Tx) {
  # check design formula: does Tx exist, is it a factor
  # display a level and control vs. contrasts

    res <- lapply(object, .wrap_limma, design)
    names(res) <- names(object)

    return(res)
  }
)

#' @rdname talus_limma
#' @aliases talus_limma,TalusDataSet
#' @exportMethod talus_limma
setMethod("talus_limma", signature(object = "TalusDataSet"),
  function(object, design = ~ 0 + Tx) {
  # check design formula: does Tx exist? is it a factor?
  # display a level and control vs. contrasts
    res <- .wrap_limma(object, design)
    return(res)
  }
)


.wrap_limma <- function(object, design) {
  # assemble row_data to be join to the output of topTable()
  row_data <- as.data.frame(rowData(object)) %>%
    tibble::rownames_to_column(var = "id")

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
  contrast_matrix <- do.call(
    makeContrasts,
    c(
      as.list(contrast_exprs),
      list(levels = mm)
    )
  )

  # 7. Fit contrasts & compute moderated statistics
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  # make data.frame
  top_table <- lapply(names(contrast_exprs), function(conts) {
    topTable(fit2,
      coef = conts, number = nrow(assay),
      sort.by = "P"
    ) %>%
      tibble::rownames_to_column(var = "id") %>%
      dplyr::left_join(row_data, by = "id")
  })
  names(top_table) <- names(contrast_exprs)

  return(top_table)
}
