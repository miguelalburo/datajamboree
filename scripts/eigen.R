#
# Eigengene X Eigenmetabolite Correlations
# 
#


# Setup -------------------------------------------------------------------

# Essentials
source("functions.R")

# Specifics
library(WGCNA)


# Load data

## Eigengenes
egenes <- read.csv(file = "../results/wgcna/eigengenes.csv", row.names = 1)
colnames(egenes) <- gsub('ME', 'G_', colnames(egenes))

## Eigenmetabolites
emets <- read.csv(file = "../results/wmcna/eigenmetabolites.csv", row.names = 1)
colnames(emets) <- gsub('ME', 'M_', colnames(emets))

## Shared Samples Filter
shared_ids <- intersect(rownames(egenes), rownames(emets))
emets <- emets[shared_ids,]
egenes <- egenes[shared_ids,]


# Correlation -------------------------------------------------------------

# Correlation Matrix
module_trait_cor  <- cor(egenes, emets, use = "p")

# P values
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(egenes))
module_trait_padj <- matrix(
  p.adjust(as.vector(module_trait_pval), method = "BH"),
  nrow      = nrow(module_trait_pval),
  dimnames  = dimnames(module_trait_pval)
)

# Label significant module-trait
text_matrix <- ifelse(
  abs(module_trait_cor) >= 0.25 & module_trait_padj <= 0.05,
  paste0(signif(module_trait_cor, 2), "\n(", signif(module_trait_padj, 1), ")"),
  NA
)
dim(text_matrix) <- dim(module_trait_cor)


# Plotting Heatmap --------------------------------------------------------


pdf(file = "../results/eigen/eigen_heatmap.pdf", width = 12, height = 8)
par(mar = c(6, 9, 3, 3))

labeledHeatmap(
  Matrix        = module_trait_cor,
  xLabels       = colnames(module_trait_cor),
  yLabels       = rownames(module_trait_cor),
  ySymbols      = rownames(module_trait_cor),
  xSymbols      = colnames(module_trait_cor),
  colorLabels   = FALSE,
  colors        = blueWhiteRed(50),
  textMatrix    = text_matrix,
  setStdMargins = FALSE,
  cex.text      = 0.5,
  zlim          = c(-1, 1),
  main          = "Module-Trait Relationships (r, BH-adj. p-value)"
)

dev.off()
