# talusR — Lightweight R toolkit for Talus Bio proteomics

For a full walk-through, please veiw the [vignette](https://chaochaowong.github.io/talusR/).

## Overview
__talusR__ provides a streamlined workflow for Talus mass-spectrometry data: import, QC, differential analysis with _limma_, and visualization. It requires DIA-NN outputs and matching sample metadata.

## Key Features
- **Import & Structure** 
  - Load DIA-NN outputs into S4 containers: TalusDataSet or TalusDataSetList for __cytoplasmic__, __nucleoplasmic__, and __chromatin__-bound fractions .
  - Built on Bioconductor classes for interoperability.
- **QC** 
  PCA and multivariate checks to assess clustering and detect outliers.
- **Differential Analysis**  
  Thin wrappers around `limma` and `matrixTests` for hypothesis testing.
- **Visualization**  
  PCA, volcano, and per-protein plots powered by `ggplot2`.
- **Bioconductor Foundations**  
  Uses `SummarizedExperiment` and `S4Vectors` under the hood.