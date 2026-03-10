# ── CRAN packages ────────────────────────────────────────────────
cran_packages <- c(
  "igraph",
  "dynamicTreeCut",
  "fastcluster",
  "WGCNA",
  "reshape2",
  "tidyverse",
  "dplyr",
  "stringr",
  "matrixStats",
  "shiny",
  "umap",
  "Rtsne",
  "caret",
  "mlr3",
  "kernlab",
  "randomForest",
  "factoextra",
  "pheatmap",
  "gplots",
  "ggraph",
  "plotly",
  "VennDiagram",
  "iml",
  "lime",
  "gprofiler2"
)

install.packages(cran_packages, repos = "https://cloud.r-project.org")

# ── Bioconductor packages ─────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_packages <- c(
  "DESeq2",
  "limma",
  "sva",
  "mixOmics",
  "GSEABase",
  "omicade4",
  "MOFA2",
  "MOFAdata",
  "cqn"
)

BiocManager::install(bioc_packages)

# ── GitHub packages ───────────────────────────────────────────────
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

remotes::install_github("igordot/shapr")         # shapr
remotes::install_github("sumbose/iRF")            # iRF
remotes::install_github("AGS-/RGCCA")             # RGCCA (if not on CRAN)

message("All packages installed successfully.")