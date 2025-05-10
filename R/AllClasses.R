#' @rdname TalusDataSet
#' @docType class
#' @aliases TalusDataSet-class
#' @export
#'
# 1) Define class
setClass(
  "TalusDataSet",                    # your new class name
  contains = "SummarizedExperiment", # inherit all the SE machinery
  slots    = list(
    intensity_group = "character",    # protein or peptide
    metric = "character"             # DIA-NN
  )
)

# 2) validity
setValidity("TalusDataSet", function(object) {

  if (! ("abund" %in% assayNames(object)) )
    return( "the assays slot must contain a matrix named 'abound'" )

  ig <- object@intensity_group
  # require exactly one value, and it must be "protein" or "peptide"
  if (length(ig) != 1 || ! ig %in% c("protein","peptide")) {
    return("`intensity_group` must be a single string: either 'protein' or 'peptide'")
  }

  TRUE
})

#' TalusDataSet object and constructors
#'
#' \code{TalusDataSet} is a subclass of \code{SummerizedExperiment} used
#' to store the input protein abundance values and addition annotation
#' of intensity group (protein or peptide) and input format (DIA-NN)
#' analysis of differential expression.
#' @rdname TalusDataSet
#' @export
#'
TalusDataSet <- function(abund_data, col_data, row_data,
                         intensity_group = 'protein',
                         metric = 'DIA-NN')

{
  # check that these agree in number
  stopifnot(ncol(abund_data) == nrow(col_data))
  stopifnot(nrow(abund_data) == nrow(row_data))

  # we expect a matrix of counts, which are non-negative integers
  abund_ata <- as.matrix(abund_ata )

  if (is(col_ata, "data.frame"))
    col_data <- as(col_data, "DataFrame")

  if (is(row_ata, "data.frame"))
    row_data <- as(row_data, "DataFrame")

  # check if the rownames of colData are simply in different order
  # than the colnames of the countData, if so throw an error
  # as the user probably should investigate what's wrong
  if (!is.null(rownames(col_data)) & !is.null(colnames(abund_data))) {
    if (all(sort(rownames(col_data)) == sort(colnames(abound_data)))) {
      if (!all(rownames(col_data) == colnames(abound_data))) {
        stop(paste("rownames of the col_data:
  ",paste(rownames(col_data), collapse=","), "
  are not in the same order as the colnames of the countData:
  ",paste(colnames(abund_ata), collapse=",")))
      }
    }
  }

  if (is.null(rownames(col_data)) & !is.null(colnames(abund_data))) {
    rownames(col_data) <- colnames(abund_ata)
  }

  se <- SummarizedExperiment(assays = SimpleList(abund=abund_data),
                             colData = col_data,
                             rowData = row_data)
  new("TalusDataSet", se,
      intensity_group = intensity,
      metric = metric)
}

