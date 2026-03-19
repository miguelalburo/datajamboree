#
# Weighted Metabolite Co-Abundance Network Analysis (WMCNA)
#
# - Identifies modules of co-abundant metabolites
# - Correlates modules with experimental traits and chemical concentrations
# - Exports a combined edge list for significant modules
#


# Setup -------------------------------------------------------------------

source("functions.R")

library(WGCNA)
library(igraph)
library(dplyr)
library(tidyr)
library(tibble)

options(stringsAsFactors = FALSE)
# enableWGCNAThreads()  # comment out in RStudio

set.seed(67)


# Parameters --------------------------------------------------------------

params <- list(
  
  # Soft thresholding
  powers           = 2:20,
  r2_cutoff        = 0.80,      # relaxed from 0.90: scale-free topology is
  # less strict for metabolite co-abundance networks
  soft_power       = NULL,      # set manually to override auto-selection
  
  # Module detection
  min_module_size  = 20,
  deep_split       = 2,
  merge_cut_height = 0.25,
  
  # Module-trait filtering
  cor_threshold    = 0.25,
  p_threshold      = 0.05,
  
  # Edge export: minimum TOM weight to include an edge
  tom_threshold    = 0.25,
  
  # Relative output directory path
  out_dir          = "../results/wmcna"
)

dir.create(params$out_dir, showWarnings = FALSE, recursive = TRUE)


# Load Data ---------------------------------------------------------------

# Metabolomics data (positive mode, glog-transformed PQN)
dims <- load.dataset(
  meta.file  = "../data/sample_sheet.csv",          meta.sep  = ",",
  data.file  = "../data/polar_pos_pqn_imputed_glog.csv", data.sep = ","
)

data      <- dims$data.matrix   # samples x metabolites
meta_data <- dims$meta.data

cat(sprintf("Loaded: %d samples x %d metabolites\n", nrow(data), ncol(data)))


# Chemical Data
chem <- read.csv(file = "../data/chem_conc.csv", row.names = 1)
cois <- read.table(file = "../results/rgcca/cois_rgcca.txt") %>% pull
# chem <- chem[, cois]

# QC: remove bad samples / metabolites ------------------------------------

gsg <- goodSamplesGenes(data, verbose = 3)
if (!gsg$allOK) {
  data <- data[gsg$goodSamples, gsg$goodGenes]
  cat(sprintf("After QC: %d samples x %d metabolites\n", nrow(data), ncol(data)))
}


# Build Trait Matrix ------------------------------------------------------

conds <- list()

# One-hot encoding site
conds$Site <- meta_data %>%
  select(Site) %>%
  rownames_to_column("rn") %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Site, values_from = value, values_fill = 0) %>%
  select(-rn)

# One-hot encoding concentration/REF
conds$REF <- meta_data %>%
  select(REF) %>%
  rownames_to_column("rn") %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = REF, values_from = value, values_fill = 0) %>%
  select(-rn)

# Chemicals
conds$Chemicals <- chem[,rownames(meta_data)] %>% t()


# Soft-Thresholding -------------------------------------------------------

sft <- pickSoftThreshold(data, powerVector = params$powers, verbose = 5)

if (is.null(params$soft_power)) {
  passed <- data.frame(sft$fitIndices) %>%
    mutate(r2 = -sign(slope) * SFT.R.sq) %>%
    filter(r2 >= params$r2_cutoff)
  
  soft_power <- if (nrow(passed) == 0) {
    warning("No power reached R2 cutoff — falling back to power 6. ",
            "Consider lowering r2_cutoff or setting soft_power manually.")
    6L
  } else {
    passed$Power[1]
  }
} else {
  soft_power <- params$soft_power
}

cat(sprintf("Soft power selected: %d\n", soft_power))

pdf(file = file.path(params$out_dir, "soft_threshold.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology R2",
     type = "n", main = "Scale Independence")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = params$powers, cex = 0.9, col = "red")
abline(h = params$r2_cutoff, col = "red",  lty = 2)
abline(v = soft_power,       col = "blue", lty = 2)
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = params$powers, cex = 0.9, col = "red")
abline(v = soft_power, col = "blue", lty = 2)
dev.off()


# Topological Overlap Matrix ----------------------------------------------

# NOTE: computationally expensive — run on HPC, then reload with readRDS()
adjacency <- adjacency(data, power = soft_power)
TOM       <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM


# Module Detection --------------------------------------------------------

metabolite_tree <- hclust(as.dist(dissTOM), method = "average")
dynamic_mods    <- cutreeDynamic(
  dendro            = metabolite_tree,
  distM             = dissTOM,
  deepSplit         = params$deep_split,
  pamRespectsDendro = FALSE,
  minClusterSize    = params$min_module_size
)
module_colors <- labels2colors(dynamic_mods)

# Merge Similar Modules ---------------------------------------------------

merge_result  <- mergeCloseModules(data, module_colors,
                                   cutHeight = params$merge_cut_height,
                                   verbose   = 3)
module_colors <- merge_result$colors
MEs           <- merge_result$newMEs

cat(sprintf("Modules after merging: %d\n", length(unique(module_colors)) - 1))
print(sort(table(module_colors), decreasing = TRUE))

pdf(file = file.path(params$out_dir, "metabolite_dendrogram.pdf"), width = 14, height = 9)
plotDendroAndColors(
  metabolite_tree,
  cbind(labels2colors(dynamic_mods), module_colors),
  c("Before merge", "After merge"),
  dendroLabels = FALSE, hang = 0.03,
  addGuide     = TRUE,  guideHang = 0.05,
  main = "Metabolite Dendrogram and Module Colours"
)
dev.off()


# Module-Trait Correlation ------------------------------------------------

# Ordered Eigenmetabolites
MEs_ordered <- orderMEs(MEs)

# Results list
full_results <- data.frame(
  module = character(),
  trait = character(),
  correlation = numeric(),
  padj = numeric()
)

# Module-trait plotting loop
for (name in names(conds)) {
  cond <- conds[[name]]
  module_trait_cor  <- cor(MEs_ordered, cond, use = "p")
  module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(data))
  module_trait_padj <- matrix(
    p.adjust(as.vector(module_trait_pval), method = "BH"),
    nrow      = nrow(module_trait_pval),
    dimnames  = dimnames(module_trait_pval)
  )
  
  # Results Table
  results <- left_join(
    as.data.frame(module_trait_cor) %>%
      rownames_to_column("module") %>%
      pivot_longer(-module, names_to = "trait", values_to = "correlation"),
    as.data.frame(module_trait_padj) %>%
      rownames_to_column("module") %>%
      pivot_longer(-module, names_to = "trait", values_to = "padj"),
    by = c("module", "trait")
  )
  
  full_results <- rbind(full_results, results)
  
  # Label significant module-trait
  text_matrix <- ifelse(
    abs(module_trait_cor) >= 0.25 & module_trait_padj <= 0.05,
    paste0(signif(module_trait_cor, 2), "\n(", signif(module_trait_padj, 1), ")"),
    NA
  )
  dim(text_matrix) <- dim(module_trait_cor)
  
  # Plotting heatmap
  pdf(file = file.path(params$out_dir, paste0(name,"_heatmap.pdf")), width = 12, height = 8)
  par(mar = c(6, 9, 3, 3))
  
  labeledHeatmap(
    Matrix        = module_trait_cor,
    xLabels       = colnames(cond),
    yLabels       = rownames(module_trait_cor),
    ySymbols      = rownames(module_trait_cor),
    colorLabels   = FALSE,
    colors        = blueWhiteRed(50),
    textMatrix    = text_matrix,
    setStdMargins = FALSE,
    cex.text      = 0.5,
    zlim          = c(-1, 1),
    main          = "Module-Trait Relationships (r, BH-adj. p-value)"
  )
  
  dev.off()
}

# Identify Significant Modules --------------------------------------------

sig_modules <- full_results %>% 
  filter(abs(correlation) > params$cor_threshold, padj < params$p_threshold)

cat(sprintf("\nSignificant modules (n = %d): %s\n",
            length(unique(sig_modules$module)), paste(unique(sig_modules$module), collapse = ", ")))

if (length(sig_modules) == 0) {
  stop("No modules passed thresholds. Check: (1) trait matrix is numeric, ",
       "(2) soft threshold plot, (3) consider relaxing cor_threshold or p_threshold.")
} else {
  write.csv(sig_modules,
            file.path(params$out_dir, "sig_modules.csv"),
            row.names = FALSE)
}
# Build Combined Edge List  -------------------------

probes    <- colnames(data)
all_edges <- list()

for (mod in unique(sig_modules$module)) {
  
  cat(sprintf("\nProcessing module: %s\n", mod))
  
  in_module  <- module_colors == gsub('ME', '', mod)
  mod_probes <- probes[in_module]
  
  # Recover TOM for this module from dissTOM
  mod_TOM           <- (1 - dissTOM)[in_module, in_module]
  dimnames(mod_TOM) <- list(mod_probes, mod_probes)
  
  # Upper triangle only (undirected, no duplicate edges)
  upper_idx <- which(upper.tri(mod_TOM) & mod_TOM >= params$tom_threshold,
                     arr.ind = TRUE)
  
  if (nrow(upper_idx) == 0) {
    warning(sprintf("Module %s: no edges above TOM threshold %.2f — skipping. ",
                    mod, params$tom_threshold),
            "Consider lowering tom_threshold.")
    next
  }
  
  mod_edges <- data.frame(
    source_metabolite = mod_probes[upper_idx[, 1]],
    target_metabolite = mod_probes[upper_idx[, 2]],
    tom_weight        = mod_TOM[upper_idx],
    module            = mod,
    stringsAsFactors  = FALSE
  )
  all_edges[[mod]] <- mod_edges
  
  cat(sprintf("  Metabolites: %d | Edges: %d=\n",
              length(mod_probes), nrow(mod_edges)))
}


# Save Outputs ------------------------------------------------------------


# Edge List
edge_list <- bind_rows(all_edges)
write.csv(edge_list,
          file.path(params$out_dir, "combined_edge_list.csv"),
          row.names = FALSE)

# Module eigenmetabolites
write.csv(MEs_ordered, row.names = T,
          file.path(params$out_dir, "eigenmetabolites.csv"))


# Save Metabolite Module lists
source("annotation.R")
for (module in names(all_edges)) {
  metabolites <- unique(c(all_edges[[module]]$source_metabolite,all_edges[[module]]$target_metabolite))
  
  metabs <- mz_to_metabolite(metabolites, pos.kegg)
  
  file_name <- paste0("../results/wmcna/annotated_metabolite_lists/",module,"_KEGG.txt")
  
  write.table(metabs, file_name,
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
}


# Summary -----------------------------------------------------------------

cat("\n--- Output Summary ---\n")

cat(sprintf("combined_edge_list.csv : %d edges across %d modules\n",
            nrow(edge_list), length(all_edges)))

cat(sprintf("Output directory       : %s\n", params$out_dir))

cat("\nEdges per module:\n")
edge_list %>% group_by(module) %>% count(module, name = "n_edges") %>% print()
