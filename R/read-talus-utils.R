# .make_se() and .construct_se()
.make_se <- function(tb,
                     meta,
                     split_by_fraction,
                     which_run,
                     which_fraction,
                     which_proteinid,
                     intensity_group = "protein",
                     metric = "DIA-NN") {
  all_runs <- meta %>%
    pull(which_run)

  if (split_by_fraction) {
    se <- meta %>%
      split(.[[which_fraction]]) %>%       # named list by Fraction values
      map(function(sub_meta) {
        .construct_tds(sub_meta,
                      which_run,
                      tb,
                      which_proteinid,
                      all_runs)
      })

    # TODO? clean up rownames: one single protein.ID ??

  }


  if (!split_by_fraction) {
    se <- .construct_tds(meta,
                        which_run,
                        tb,
                        which_proteinid,
                        all_runs,
                        intensity_group,
                        metric)

  }

  return(se)
}



.construct_tds <- function(sub_meta,
                          which_run,
                          tb,
                          which_proteinid,
                          all_runs,
                          intensity_group = "protein",
                          metric = "DIA-NN") {
  offset = 1
  runs <- sub_meta[[which_run]]
  # assay
  assay <- tb %>% dplyr::select(any_of(runs)) %>%
    as.matrix()
  rownames(assay) <- tb %>% pull(which_proteinid)

  # colData
  col_data <- as(sub_meta, 'DataFrame')
  rownames(col_data) <- runs

  # rowData
  row_data <- tb %>% dplyr::select(-any_of(all_runs))
  row_data <- as(row_data, 'DataFrame')
  rownames(row_data) <- tb %>% pull(which_proteinid)

  se <- TalusDataSet(assay_data = assay,
                     col_data = col_data,
                     row_data = row_data,
                     intensity_group = intensity_group,
                     metric = metric)

  #se <- SummarizedExperiment(assays=log2(assay+offset),
  #                           rowData=row_data,
  #                           colData=col_data,
  #                           metadata=col_data,
  #                           checkDimnames=TRUE)
}

.remove_few_measurements <- function(x, threshold = 0.85) {
  na_frac <- rowMeans(is.na(assays(x)[[1]]))
  x <- x[na_frac < threshold]
}
