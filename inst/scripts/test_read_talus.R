# test_read_talus.R
library(SummarizedExperiment)
library(rlang)
library(dplyr)
library(purrr)
library(ggplot2)

# devtools::install_github('chaochaowong/talusR')
# library(talusR)

work_dir <- file.path('///Volumes/sarthy_j/Sarthy/Sarthy_Lab/Proteomics')
thp1_dir <- file.path(work_dir, 'talus_results/GB_THP1_anthras_20240716')
data_dir <- file.path(work_dir, 'talus/GB_THP1_anthras_20240716/data')
stat_dir <- file.path(thp1_dir, 'stats')
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

# split by Frx
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

plot_pca(se_frac, top_n = 1000, color_by = 'Tx')

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
res_t <- talus_row_t_welch(se_frac, design = ~ 0 + Tx)

# summary_talus() can only do one fraction at once
chrom_filt <-
  summary_talus(res_t[['chrom']],
                alpha = 0.05,
                lfc_threshold = 0.5)

writexl::write_xlsx(res_t[['chrom']],
                    path = here::here(stat_dir,
                                      'chrom-treatment-vs-DMSO-ttest.xlsx'))
writexl::write_xlsx(chrom_filt,
                    path = here::here(stat_dir,
                                      'chrom-treatment-vs-DMSO-ttest-significance-lfc0.5-alpha5e-2.xlsx'))
#
#
# limma results
#
res_limma <- talus_limma(se_frac, design = ~ 0 + Tx)
chrom_filt <-
  summary_talus(res_limma[['chrom']],
                alpha = 0.05,
                lfc_threshold = 0.5)

writexl::write_xlsx(res_limma[['chrom']],
                    path = here::here(stat_dir,
                                      'chrom-treatment-vs-DMSO-limma.xlsx'))
writexl::write_xlsx(chrom_filt,
                    path = here::here(stat_dir,
                                      'chrom-treatment-vs-DMSO-limma-significance-lfc0.5-alpha5e-2.xlsx'))
#
# per protein plot
#
protein_id = 'Q969X6;Q969X6-2;Q969X6-3'
#protein_id = "A0AVT1;A0AVT1-3;A0AVT1-4"
per_protein_abun(se_frac, protein_id = protein_id,
                 category_by = 'Tx')

per_protein_abun(se_frac[['chrom']], protein_id = protein_id,
                 category_by = 'Tx') +
  labs(title='UTP4')

#
# plot_volcano?
#

#
# differential abundance vs. DepMap
#
