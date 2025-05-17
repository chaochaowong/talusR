#' Import Talus files
#'
#' Import Talus files that were processed by DIANN and convert to
#' SummarizedExperiment instance. Can split the data into three instancs
#' containing fractions in nuc, cyto, and plasm.
#'
#' @param file a path to a mass-spectrum intensity signal matric file (TSV OR CSV) or a connection.
#' @param meta_data a path to the metadata file (CSV) for the runs (samples).
#' @param which_proteinid which column in the \code{file} input represents the protein id
#' @param which_sequence which column in the \code{file} input represents the sequence of the protein
#' @param which_run which column in the meta file represents the runs (sample) ID.
#' @param remove_few_measurements if \code{TRUE}, remove measurements that are absent 85% across all samples
#' @param split_by_fraction if \code{TRUE}, split the protein abundance file according to value indicating in the values of the \code{which_fraction} column in the meta file
#' @param log_transform if \code{TRUE}, perform log2 transformation to the protein abundance
#' @param intensity_group  defual to "protein".
#' @param metric default to "DIA-NN".
#' @return a \code{SummarizedExperiment} objects with three layers of assay data
#'
#' @importFrom rlang has_name
#' @importFrom readr read_delim read_csv
#' @importForm tools file_ext
#' @importFrom purrr map
#' @import SummarizedExperiment
#' @import dplyr
#' @examples
#'
#' @export
read_talus <- function(file, meta_file,
                       which_proteinid = "Protein.Ids",
                       which_fraction = "Frx",
                       which_sequence = NA,
                       which_run = "Run",
                       remove_few_measurements = TRUE,
                       split_by_fraction = TRUE,
                       log_transform = TRUE,
                       rowname_repair = TRUE,
                       intensity_group = "protein",
                       metric = "DIA-NN") {
  # assume 1) file is in tsv format 2) meta_data is in csv format
  require(readr)
  require(SummarizedExperiment)
  require(dplyr)

  #
  # 1. read intensity signal
  #
  is_csv  <- tolower(tools::file_ext(file)) == "csv"
  is_tsv  <- tolower(tools::file_ext(file)) == "tsv"

  if (!any(is_csv, is_tsv))
    stop('File must be in either csv or tsv format.')

  if (is_tsv) {
    tb <- read_delim(file, delim = "\t")
  }

  if (is_csv) {
    tb <- read_csv(file)
  }

  #
  # 2. read meta_file
  #
  if (tolower(tools::file_ext(meta_file)) == "csv")
    meta <- read_csv(meta_file)
  else
    stop('Meta file must be in csv format.')

  #
  # 3. check which_run column exists in meta
  #
  if (!rlang::has_name(meta, which_run)) {
    stop(sprintf("which_run: column '%s' not found in %s", which_run, meta_file))
  }

  #
  # 4. check which_proteinid column exists
  #
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
  tds <- .make_tds(
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
