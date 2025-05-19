# .make_se() and .construct_se()
.make_tds <- function(tb,
                      meta,
                      split_by_fraction,
                      which_run,
                      which_fraction,
                      which_proteinid,
                      intensity_group = "protein",
                      metric = "DIA-NN",
                      log_transform = 'log2') {

  all_runs <- meta %>%
    pull(which_run)

  if (split_by_fraction) {
    tds <- meta %>%
      split(.[[which_fraction]]) %>%       # named list by Fraction values
      map(function(sub_meta) {
        .construct_tds(sub_meta,
                       which_run,
                       tb,
                       which_proteinid,
                       all_runs,
                       intensity_group,
                       metric,
                       log_transform)
      })

    tds <- TalusDataSetList(tds)
    # TODO? clean up rownames: one single protein.ID ??

  }


  if (!split_by_fraction) {
    tds <- .construct_tds(meta,
                          which_run,
                          tb,
                          which_proteinid,
                          all_runs,
                          intensity_group,
                          metric,
                          log_transform)

  }

  return(tds)
}

.construct_tds <- function(sub_meta,
                           which_run,
                           tb,
                           which_proteinid,
                           all_runs,
                           intensity_group = "protein",
                           metric = "DIA-NN",
                           log_transform = 'log2') {
  offset = 1
  runs <- sub_meta[[which_run]]
  # assay
  assay <- tb %>% dplyr::select(any_of(runs)) %>%
    as.matrix()
  assay <- log2(assay + offset)
  rownames(assay) <- tb %>% pull(which_proteinid)

  # colData
  col_data <- as(sub_meta, 'DataFrame')
  rownames(col_data) <- runs

  # rowData
  row_data <- tb %>% dplyr::select(-any_of(all_runs))
  row_data <- as(row_data, 'DataFrame')
  rownames(row_data) <- tb %>% pull(which_proteinid)

  TalusDataSet(assay_data = assay,
               col_data = col_data,
               row_data = row_data,
               intensity_group = intensity_group,
               metric = metric,
               log_transform = log_transform)
}

.remove_few_measurements <- function(x, threshold = 0.85) {
  na_frac <- rowMeans(is.na(assays(x)[[1]]))
  x <- x[na_frac < threshold]
}
