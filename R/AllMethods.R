setMethod(
  "show", "TalusDataSet",
  function(object) {
    ## First, show the SummarizedExperiment summary:
    callNextMethod()

    ## Then fetch and print your metadata fields:
    md <- metadata(object)
    ig <- if ("intensity_group" %in% names(md)) md$intensity_group else NA
    mt <- if ("metric"          %in% names(md)) md$metric          else NA
    lt <- if ("log_transform"   %in% names(md)) md$log_transform   else NA


    cat("intensity_group:", ig, "\n")
    cat("metric:         ", mt, "\n")
    cat("log_transform:  ", lt, "\n")
  }
)

#' @rdname TalusDataSetList
#' @importMethodsFrom S4Vectors show
setMethod("show", "TalusDataSetList", function(object) {
  ## 1) show the SimpleList summary:
  callNextMethod()

  ## 2) then show each TalusDataSet inside:
  nms <- names(object)
  for (i in seq_along(object)) {
    header <- if (!is.null(nms) && nzchar(nms[i])) {
      paste0("[[", i, "]] — name: ", nms[i], "\n")
    } else {
      paste0("[[", i, "]]\n")
    }
    cat("\n", header, sep = "")
    show(object[[i]])
  }
})
