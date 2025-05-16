# tmp file just for testing
work_dir <- '///Volumes/sarthy_j/Sarthy/Sarthy_Lab/Proteomics'
thp1_dir <- here::here(work_dir,
                       'talus_results/GB_THP1_anthras_20240716')
data_dir <- here::here(work_dir,
                       'talus/GB_THP1_anthras_20240716/data')
stat_dir <- here::here(thp1_dir, 'stats')
fig_dir  <- here::here(thp1_dir, 'figures')
# file path
protein_file <- here::here(data_dir,
                           'THP1_proteinIDs_matrix.tsv')
meta_file    <- here::here(data_dir,
                           'THP1_meta_data.csv')
file = protein_file
meta = meta_file
which_proteinid = 'Protein.Ids'
which_fraction = 'Frx'
which_sequence = NA
which_run = 'Run'
remove_few_measurements = TRUE
split_by_fraction = TRUE


# for each fraction; select 500 protein
tb <- read_delim(file, delim = "\t")

# only keep 'DMSO, Doxo, and Etopo'
meta <- read_csv(meta_file)

#
# slice proteins and samples
#
meta_sub <- meta %>%
  dplyr::filter(Tx %in% c('DMSO', 'Doxo', 'Etopo'))

set.seed(42)  # for reproducibility
tb_sub <- tb %>%
  slice_sample(n = 1500) %>%
  dplyr::select(names(tb)[1:5],
                meta_sub$Run)

#
# make the following test data sets:
# - extdata/test_data.tsv
# - extdata/test_meta.csv
#  - data/test_se.rda
#





