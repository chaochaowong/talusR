# AllMethods.R()
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
