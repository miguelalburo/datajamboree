#
# Regularised Generalised Canonical Correlation Analysis (RGCCA)
#
# Chemicals ~ RNAseq Gene Counts ~ DIMS Peak Intensity
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


# Naming blocks for clarity in outputs
block.names <- c("Chemicals", "Transcriptomics", "Metabolomics")

# Scaling
A <- lapply(A, function(x) scale(x))


# Design Matrix
# Chemicals -> Transcriptomics -> Metabolomics
# Chemicals and Metabolomics are NOT directly connected
C <- matrix(c(0, 1, 1,
              1, 0, 1,
              1, 1, 0), 3, 3)


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


# Features-of-Interest ----------------------------------------------------


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
      
      name <- paste(cond, p1, p2, sep = "_")
      
      p[[name]] <- factor.plot(y1, y2, group) +
        labs(title = paste(p1, "vs", p2),
             color = cond,
             x = p1,
             y = p2)
    }
  }
}

library(plotly)

conds <- c("Site", "REF", "Description")
p.3D <- list()

for (cond in conds) {
  for (component in seq(ncomp)) {
    
    group <- rna.ds$meta.data[shared.ids, ][[cond]]
    
    # One axis per block
    y1 <- cca.res$Y[[1]][, component]  # Chemicals
    y2 <- cca.res$Y[[2]][, component]  # Transcriptomics
    y3 <- cca.res$Y[[3]][, component]  # Metabolomics
    
    name <- paste(cond, "component", component, sep = "_")
    
    p.3D[[name]] <- plot_ly(
      x = y1, y = y2, z = y3,
      color = group,
      colors = "Set1",
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 5, opacity = 0.8),
      text = paste("Sample:", shared.ids,
                   "<br>Group:", group),
      hoverinfo = "text"
    ) %>%
      layout(
        title = paste(cond, "- Component", component),
        scene = list(
          xaxis = list(title = "Chemicals"),
          yaxis = list(title = "Transcriptomics"),
          zaxis = list(title = "Metabolomics")
        )
      )
    
    # Save 3D plots
    file = paste0("../results/rgcca/factor_plots_3D_", name, ".html")
    htmlwidgets::saveWidget(p.3D[[name]], file)
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
  pdf(paste('../results/rgcca/factor_plots_', fplot, '.pdf', sep = ''),
      width = 10, height = 6)
  print(p[[fplot]])
  dev.off()
}


# Write Block Weights to csv
for (i in seq_along(cca.res$a)) {
  block.name <- block.names[i]
  weights <- cca.res$a[[i]] %>% 
    as.data.frame() %>% filter(V1 != 0, V2 != 0)
  out.file <- paste0('../results/rgcca/block_weights_', block.name, '.csv')
  write.csv(weights, out.file, row.names = TRUE)
}


# Write to report

out.path <- '../results/rgcca/rgcca_report.txt'
sink(out.path)

cat("════════════════════════════════════════════════════════════\n")
cat("  RGCCA Results\n")
cat("════════════════════════════════════════════════════════════\n\n")

## Parameters
cat("── Model Parameters ────────────────────────────────────────\n")
cat(sprintf("  Number of components : %d\n", ncomp))
cat(sprintf("  Tau                  : %s\n", paste(tau, collapse = ", ")))
cat(sprintf("  Sparsity             : %s\n", paste(sparsity, collapse = ", ")))
cat("\n  Connection Matrix (C):\n")
C.labelled <- C
rownames(C.labelled) <- colnames(C.labelled) <- block.names
print(as.data.frame(C.labelled))
cat("\n")


## AVE Outer (Variance Explained per Block) 
cat("── AVE Outer: Variance Explained per Block ─────────────────\n")
ave.outer <- cca.res$AVE$AVE_X  # Named list, one entry per block
for (i in seq_along(ave.outer)) {
  cat(sprintf("\n  Block %d | %s\n", i, block.names[i]))
  for (comp in seq_along(ave.outer[[i]])) {
    cat(sprintf("    Component %d : %.4f (%.2f%%)\n",
                comp,
                ave.outer[[i]][comp],
                ave.outer[[i]][comp] * 100))
  }
}
cat("\n")


## AVE Inner
cat("── AVE Inner: Inter-block Variance Explained ───────────────\n")
ave.inner <- cca.res$AVE$AVE_inner  # One value per component
for (comp in seq_along(ave.inner)) {
  cat(sprintf("  Component %d : %.4f (%.2f%%)\n",
              comp,
              ave.inner[comp],
              ave.inner[comp] * 100))
}
cat("\n")

cat("════════════════════════════════════════════════════════════\n")
cat(sprintf("  File written: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("════════════════════════════════════════════════════════════\n")

sink()
message("Results written to: ", out.path)



