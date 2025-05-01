#' Import Talus files
#'
#' Import Talus files that were processed by Diann and convert to
#' SummarizedExperiment instances with three layers: nuc, cyto, and plasm.
#'
#' @param file A path to a Diann summarized protein file or a connection.
#' @param meta_data A path to a meta (column) data file for samples.
#' @param which_proteinid
#' @param which_sequence
#' @param which_run
#' @param remove_few_measurements If \code{TRUE}, remove measurements that are absent 85% across all samples
#' @param split_by_fraction If \code{TRUE}, split the protein abundance file according to value indicating in the values of the \code{which_fraction} column in the meta file
#' @param log_transform If \code{TRUE}, perform log2 transformation to the protein abundance
#'
#' @return A \code{SummarizedExperiment} objects with three layers of assay data
#' @export
read_talus <- function(file, meta_file,
                       which_proteinid = 'Protein.Ids',
                       which_fraction = 'Frx',
                       which_sequence = NA,
                       which_run = 'Run',
                       remove_few_measurements = TRUE,
                       split_by_fraction = TRUE,
                       log_transform = TRUE,
                       rowname_repair = TRUE
                       ) {

  # assume 1) file is in tsv format 2) meta_data is in csv format
  require(readr)
  require(SummarizedExperiment)

  tb <- read_delim(file, delim='\t')
  meta <- read_csv(metadata)

  # check if the column defined by which_run exist in meta
  if (! rlang::has_name(meta, which_run)) {
    stop(sprintf("Column '%s' not found in %s", which_run, meta))
  }

  if (!rlang::has_name(tb, which_proteinid)) {
    stop(sprintf("Column '%s' not found in %s", which_proteinid, file))
  }


  # get run ID
  runs <- meta %>%
    pull(which_run)

  # check
  if (!all(rlang::has_name(tb, runs))) {
    stop(sprintf("Runs ID %s not found in %s",
                 paste(runs[!has_name(tb, runs)], collapse = ','),
                 file))
  }


  #
  # convert protein abundance to a list of SumarizedExperiment instance
  #
  se <- .make_se(tb, meta, split_by_fraction,
                 which_run, which_fraction)

  #
  # remove_few_measurements
  #
  if (remove_few_measurements) {
    if (is.list(se))
      se <- lapply(se, .remove_few_measurements, threshold = 0.85)
    else
      se <- .remove_few_measurements(se, threshold = 0.85)
  }

  #
  # tidy up rownames
  #

  return(se)
}
