## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
opts_chunk$set(
  eval = TRUE, 
  echo = TRUE, 
  warning = FALSE, 
  message = FALSE,
  fig.align = 'center',
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)

## ----install-talusR, eval=FALSE-----------------------------------------------
# devtools::install_github('chaochaowong/talusR',
#                          ref          = "dev",
#                          dependencies = TRUE,
#                          build_vignettes = TRUE
# )
# 
# # or to install the Bioconductor dependencies
# BiocManager::install(
#   "chaochaowong/talusR@dev",
#   dependencies = TRUE
#   build_vignettes = TRUE
# )

## ----load-library-------------------------------------------------------------
library(talusR)

## ----assemble-files-----------------------------------------------------------
# get paths
protein_file <- system.file('extdata', 'test_data.tsv', 
                            package = 'talusR')
meta_file    <- system.file('extdata', 'test_meta.csv',
                            package = 'talusR')

## ----import-test-data---------------------------------------------------------
# without split by fraction -> return TalusDataSet
tds  <- talusR::read_talus(file = protein_file,
                           meta_file = meta_file,
                           which_proteinid = "Protein.Ids",
                           which_run = "Run",
                           remove_few_measurements = TRUE,
                           split_by_fraction = FALSE)
tds

## ----pca-tds, fig.cap='PCA on TalusDataSet.'----------------------------------
talusR::plot_pca(tds, color_by='Frx')

## ----split-by-fraction--------------------------------------------------------
tdsl <- read_talus(file = protein_file,
                   meta_file = meta_file,
                   which_proteinid = "Protein.Ids",
                   which_fraction = "Frx",
                   which_sequence = NA,
                   which_run = "Run",
                   remove_few_measurements = TRUE,
                   split_by_fraction = TRUE,
                   log_transform = 'log2',
                   intensity_group = "protein",
                   metric = "DIA-NN")
# not sure why tdsl is a list!!! need to debug
tdsl

## ----pca-for-TalusDataSetList, fig.cap='PCA performed separately on each fraction in the TalusDataSetList.'----
plot_pca(tdsl, color_by='Tx')

## ----limma--------------------------------------------------------------------
res_limma <- talus_limma(tdsl, design = ~ 0 + Tx)

## ----summarize-limma----------------------------------------------------------
limma_filt <- summary_talus(res_limma[['chrom']], 
                            alpha = 0.05, 
                            lfc_threshold = 0.5)
limma_filt

## ----volcano-limma------------------------------------------------------------
plot_volcano(res_limma[['chrom']], 
             alpha = 0.05,
             lfc_threshold = 0.5, 
             use_adjP = TRUE, 
             label_top_n=10)

## ----row-t-welch, eval=FALSE--------------------------------------------------
# res_t <- talus_row_t_welch(tdsl, design = ~ 0 + Tx)
# names(res_t)

