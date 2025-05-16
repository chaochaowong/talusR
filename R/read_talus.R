#' Import Talus files
#'
#' Import Talus files that were processed by DIANN and convert to
#' SummarizedExperiment instance. Can split the data into three instancs
#' containing fractions in nuc, cyto, and plasm.
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
#'
#' @importFrom rlang has_name
#' @importFrom readr read_delim read_csv
#' @importFrom purrr map
#' @import SummarizedExperiment
#' @import dplyr
#' @export
read_talus <- function(file, meta_file,
                       which_proteinid = "Protein.Ids",
                       which_fraction = "Frx",
                       which_sequence = NA,
                       which_run = "Run",
                       remove_few_measurements = TRUE,
                       split_by_fraction = TRUE,
                       log_transform = TRUE,
                       rowname_repair = TRUE) {
  # assume 1) file is in tsv format 2) meta_data is in csv format
  require(readr)
  require(SummarizedExperiment)
  require(dplyr)

  tb <- read_delim(file, delim = "\t")
  meta <- read_csv(meta_file)

  # check which_run column exists in meta
  if (!rlang::has_name(meta, which_run)) {
    stop(sprintf("which_run: column '%s' not found in %s", which_run, meta_file))
  }

  # check which_proteinid column exists
  if (!rlang::has_name(tb, which_proteinid)) {
    stop(sprintf("which_proteinid: column '%s' not found in %s", which_proteinid, file))
  }

  # if split_by_fraction, check which_fraction column exists
  if (split_by_fraction) {
    if (!rlang::has_name(meta, which_fraction)) {
      top(sprintf("which_fraction: column '%s' not found in %s", which_fraction, meta_file))
    }
  }

  # get run ID
  runs <- meta %>%
    pull(which_run)

  # check if run in meta match what's in file (tb)
  if (!all(rlang::has_name(tb, runs))) {
    stop(sprintf(
      "Runs ID %s not found in %s",
      paste(runs[!has_name(tb, runs)], collapse = ","),
      file
    ))
  }


  #
  # convert protein abundance to a list of SumarizedExperiment instance
  #
  se <- .make_se(
    tb = tb,
    meta = meta,
    split_by_fraction = split_by_fraction,
    which_run = which_run,
    which_fraction = which_fraction,
    which_proteinid = which_proteinid
  )

  #
  # remove_few_measurements
  #
  if (remove_few_measurements) {
    if (is.list(se)) {
      se <- lapply(se, .remove_few_measurements, threshold = 0.85)
    } else {
      se <- .remove_few_measurements(se, threshold = 0.85)
    }
  }

  #
  # tidy up rownames
  #

  return(se)
}
