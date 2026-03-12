#
# Regularised Generalised Canonical Correlation Analysis (RGCCA)
#
# Chemicals -> RNAseq Gene Counts -> DIMS Peak Intensity
#


# Setup -------------------------------------------------------------------

source('functions.R')

library(RGCCA)
library(plotly)
source('fileio.R')


# Load data

## Chemicals
chem.ds <- load.dataset(
  meta.file = '../data/sample_sheet.csv', meta.sep = ',',
  data.file = '../data/chem_conc.csv', data.sep = ','
)

## Gene Counts | Normalised
rna.ds <- load.dataset(
  meta.file = '../data/sample_sheet.csv', meta.sep = ',',
  data.file = '../data/rna_norm_counts.csv', data.sep = ','
)

## DIMS Peak Intensity | Positive Ion Mode
polar.ds <- load.dataset(
  meta.file = '../data/sample_sheet.csv', meta.sep = ',',
  data.file = '../data/polar_pos_pqn_imputed_glog.csv', data.sep = ','
)


# Combining into Feature Matrix List
# Intersect across all three blocks
shared.ids <- intersect(
  intersect(row.names(chem.ds$data.matrix), row.names(rna.ds$data.matrix)),
  row.names(polar.ds$data.matrix)
)

A <- list(
  chem.ds$data.matrix[shared.ids, ],   # Block 1: Chemicals
  rna.ds$data.matrix[shared.ids, ],    # Block 2: Transcriptomics
  polar.ds$data.matrix[shared.ids, ]   # Block 3: Metabolomics
)


# Analysis ----------------------------------------------------------------


# Scaling
A <- lapply(A, function(x) scale(x))


# Design Matrix
# Chemicals -> Transcriptomics -> Metabolomics
# Chemicals and Metabolomics are NOT directly connected
C <- matrix(c(0, 1, 0,
              1, 0, 1,
              0, 1, 0), 3, 3)


# Perform RGCCA
sparsity <- c(0.2, 0.2, 0.2)   # One value per block
tau <- "optimal"
ncomp <- 2

cca.res <- rgcca(
  blocks     = A,
  connection = C,
  sparsity   = sparsity,
  tau        = tau,
  ncomp      = ncomp,
  method     = "rgcca",
  scale      = FALSE, verbose = FALSE
)


# Results -----------------------------------------------------------------


# Average Variance Explained
# Block 1 = Chemicals, Block 2 = RNA, Block 3 = Metabolomics
cca.res$AVE


# Number of features selected per block
colSums(cca.res$a[[1]] != 0)  # Chemicals
colSums(cca.res$a[[2]] != 0)  # Transcriptomics
colSums(cca.res$a[[3]] != 0)  # Metabolomics


# Variables of interest (UNION across components) per block

chemicals_of_interest <- cca.res$a[[1]] %>%
  as.data.frame() %>%
  rownames_to_column('chemical') %>%
  filter(V1 != 0 | V2 != 0) %>%
  pull(chemical)

genes_of_interest <- cca.res$a[[2]] %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  filter(V1 != 0 | V2 != 0) %>%
  pull(gene)

metabolites_of_interest <- cca.res$a[[3]] %>%
  as.data.frame() %>%
  rownames_to_column('metabolite') %>%
  filter(V1 != 0 | V2 != 0) %>%
  pull(metabolite)


# Plots -------------------------------------------------------------------


# Factor plots for all connected block pairs
# Pair A: Chemicals (1) vs Transcriptomics (2)
# Pair B: Transcriptomics (2) vs Metabolomics (3)
# Pair C: Chemicals (1) vs Metabolomics (3) — indirect, exploratory

blocks <- c("Chemicals", "Transcriptomics", "Metabolomics")
block.pairs <- list(
  c(1, 2),  # Chemicals vs Transcriptomics
  c(2, 3),  # Transcriptomics vs Metabolomics
  c(1, 3)   # Chemicals vs Metabolomics (indirect)
)

conds <- c("Site", "REF", "Description")
p <- list()

for (cond in conds) {
  for (component in seq(ncomp)) {
    for (pair in block.pairs) {
      
      group <- rna.ds$meta.data[shared.ids, ][[cond]]
      y1 <- cca.res$Y[[pair[1]]][, 1]
      y2 <- cca.res$Y[[pair[2]]][, 1]
      
      p1 <- blocks[pair[1]]
      p2 <- blocks[pair[2]]
      
      name <- paste(cond, p1, p1, sep = "_")
      
      p[[name]] <- factor.plot(y1, y2, group) +
        labs(title = paste(p1, "vs", p2),
             legend = cond,
             x = p1,
             y = p2)
    }
  }
}


# Saving & Writing --------------------------------------------------------


# Write chemicals/genes/metabolites-of-interest
write.table(chemicals_of_interest, '../results/rgcca/cois_rgcca.txt',
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(genes_of_interest, '../results/rgcca/gois_rgcca.txt',
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(metabolites_of_interest, '../results/rgcca/mois_rgcca.txt',
            row.names = FALSE, col.names = FALSE, quote = FALSE)


# Write factor plots to .pdf
for (fplot in names(p)) {
  pdf(paste('../results/rgcca/factor_plots_', fplot, '.pdf', sep = ''), width = 10, height = 6)
  print(p[[fplot]])
  dev.off()
}