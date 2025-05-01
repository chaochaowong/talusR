# test_read_talus.R
library(S4Vectors)
library(SummarizedExperiment)
library(rlang)


work_dir <- file.path('///Volumes/sarthy_j/Sarthy/Sarthy_Lab/Proteomics')
thp1_dir <- file.path(work_dir, 'talus_results/GB_THP1_anthras_20240716')
data_dir <-file.path(work_dir, 'talus/GB_THP1_anthras_20240716/data')
file <- here::here(data_dir, 'THP1_proteinIDs_matrix.tsv')
meta_file <- here::here(data_dir, 'THP1_meta_data.csv')

# test starts here
se <- read_talus(file = file,
                      meta = meta_file,
                      which_proteinid = 'Protein.Ids',
                      which_fraction = 'Frx',
                      which_sequence = NA,
                      which_run = 'Run',
                      remove_few_measurements = TRUE,
                      split_by_fraction = FALSE,
                      log_transform = TRUE,
                      rowname_repair = TRUE)

plot_pca(se, top_n = 1000, color_by = 'Frx')

se_frac <- read_talus(file = file,
                 meta = meta_file,
                 which_proteinid = 'Protein.Ids',
                 which_fraction = 'Frx',
                 which_sequence = NA,
                 which_run = 'Run',
                 remove_few_measurements = TRUE,
                 split_by_fraction = TRUE,
                 log_transform = TRUE,
                 rowname_repair = TRUE)

plot_pca(se_frac, top_n = 1000)

#
# differential analysis
#

# 1. tidy up 'Tx': make sure the fist level is the control
tx_levels <- c('DMSO', 'Doxo', 'Etopo', 'Acla', 'Dim_dox')

se_frac <- lapply(se_frac, function(x) {
  colData(x)$Tx <- factor(colData(x)$Tx, levels = tx_levels)
  x
})


#
# t Welch testing results
#
res_t <- talus_row_t_welch(se, design = ~ 0 + Tx)

#
# limma results
#
res_limma <- talus_limma(se_frac, design = ~ 0 + Tx)

#
# plot_per_protein
#

#
# plot_vocano
#

