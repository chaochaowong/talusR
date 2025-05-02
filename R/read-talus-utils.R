# .make_se() and .construct_se()
.make_se <- function(tb,
                     meta,
                     split_by_fraction,
                     which_run,
                     which_fraction,
                     which_proteinid) {
  all_runs <- meta %>%
    pull(which_run)

  if (split_by_fraction) {
    se <- meta %>%
      split(.[[which_fraction]]) %>%       # named list by Fraction values
      map(function(sub_df) {
        .construct_se(sub_df,
                      which_run,
                      tb,
                      which_proteinid,
                      all_runs)
      })
    # set assays name
    for (n in names(se)) {
      names(assays(se[[n]])) <- n
    }
    # clean up rownames: one single protein.ID

  }


  if (!split_by_fraction) {
    se <- .construct_se(meta,
                        which_run,
                        tb,
                        which_proteinid,
                        all_runs)
    names(assays(se)) <- 'all'

  }

  return(se)
}



.construct_se <- function(sub_df,
                          which_run,
                          tb,
                          which_proteinid,
                          all_runs) {
  offset = 1
  runs <- sub_df[[which_run]]
  # assay
  assay <- tb %>% dplyr::select(any_of(runs)) %>%
    as.matrix()
  rownames(assay) <- tb %>% pull(which_proteinid)
  # colData
  col_data <- as(sub_df, 'DataFrame')
  rownames(col_data) <- runs
  # rowData
  row_data <- tb %>% dplyr::select(-any_of(all_runs))
  row_data <- as(row_data, 'DataFrame')
  rownames(row_data) <- tb %>% pull(which_proteinid)

  se <- SummarizedExperiment(assays=log2(assay+offset),
                             rowData=row_data,
                             colData=col_data,
                             metadata=col_data,
                             checkDimnames=TRUE)
}

.remove_few_measurements <- function(x, threshold = 0.85) {
  na_frac <- rowMeans(is.na(assays(x)[[1]]))
  x <- x[na_frac < threshold]
}
