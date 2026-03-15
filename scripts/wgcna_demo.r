#
# Weighted Gene Co-Expression Network Analysis (WGCNA)
#
# - Identifies modules of co-expressed genes
# - Correlates modules with sample traits
# - Exports a combined edge list and hub gene table for significant modules
#


# Setup -------------------------------------------------------------------

# Essentials
source("functions.R")


# Specifics
library(WGCNA)
library(igraph)

options(stringsAsFactors = FALSE)
# enableWGCNAThreads()  # comment out in RStudio

# Seed
set.seed(67)


# Parameters --------------------------------------------------------------


params <- list(
  
  # Soft thresholding
  powers           = 2:20,
  r2_cutoff        = 0.90,
  soft_power       = NULL,
  
  # Module detection
  min_module_size  = 20,
  deep_split       = 2,
  merge_cut_height = 0.25,
  
  # Module-trait filtering
  cor_threshold    = 0.25,
  p_threshold      = 0.05,
  
  # Edge export: minimum TOM weight to include an edge
  tom_threshold    = 0.25,
  
  # Hub gene: top % by degree within each module
  hub_quantile     = 0.95,
  
  out_dir          = "../results/wgcna"
)

dir.create(params$out_dir, showWarnings = FALSE, recursive = TRUE)


# Load Data ---------------------------------------------------------------


ds <- load.dataset(
  meta.file  = "../data/sample_sheet.csv",  meta.sep  = ",",
  data.file  = "../data/rna_norm_counts.csv", data.sep = ","
)

data      <- ds$data.matrix
meta_data <- ds$meta.data

cat(sprintf("Loaded: %d samples x %d genes\n", nrow(data), ncol(data)))

conds <- lapply(meta_data, function(x) as.numeric(as.factor(x)))
conds <- as.data.frame(conds, row.names = rownames(meta_data))

gsg <- goodSamplesGenes(data, verbose = 3)
if (!gsg$allOK) {
  data <- data[gsg$goodSamples, gsg$goodGenes]
  cat(sprintf("After QC: %d samples x %d genes\n", nrow(data), ncol(data)))
}


# Soft-Thresholding -------------------------------------------------------


sft <- pickSoftThreshold(data, powerVector = params$powers, verbose = 5)

if (is.null(params$soft_power)) {
  passed <- data.frame(sft$fitIndices) %>%
    mutate(r2 = -sign(slope) * SFT.R.sq) %>%
    filter(r2 >= params$r2_cutoff)
  
  soft_power <- if (nrow(passed) == 0) {
    warning("No power reached R2 cutoff — falling back to power 6.")
    6L
  } else {
    passed$Power[1]
  }
} else {
  soft_power <- params$soft_power
}

cat(sprintf("Soft power: %d\n", soft_power))

pdf(file.path(params$out_dir, "01_soft_threshold.pdf"), width = 10, height = 5)
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


### Run on HPC, then load result here
# adjacency <- adjacency(data, power = soft_power)
# TOM <- TOMsimilarity(adjacency)
# saveRDS(TOM, file.path(params$out_dir, "TOM.rds"))

TOM     <- readRDS("../data/TOM.rds")
dissTOM <- 1 - TOM


# Module Detection --------------------------------------------------------


gene_tree    <- hclust(as.dist(dissTOM), method = "average")
dynamic_mods <- cutreeDynamic(
  dendro            = gene_tree,
  distM             = dissTOM,
  deepSplit         = params$deep_split,
  pamRespectsDendro = FALSE,
  minClusterSize    = params$min_module_size
)
module_colors <- labels2colors(dynamic_mods)


# Merge Similar Modules ---------------------------------------------------


merge_result  <- mergeCloseModules(data, module_colors,
                                   cutHeight = params$merge_cut_height,
                                   verbose = 3)
module_colors <- merge_result$colors
MEs           <- merge_result$newMEs

cat(sprintf("Modules after merging: %d\n", length(unique(module_colors)) - 1))
print(sort(table(module_colors), decreasing = TRUE))

pdf(file.path(params$out_dir, "02_gene_dendrogram.pdf"), width = 14, height = 9)
plotDendroAndColors(
  gene_tree,
  cbind(labels2colors(dynamic_mods), module_colors),
  c("Before merge", "After merge"),
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = "Gene Dendrogram and Module Colours"
)
dev.off()


# Module-Trait Correlation ------------------------------------------------


MEs_ordered       <- orderMEs(MEs)
module_trait_cor  <- cor(MEs_ordered, conds, use = "p")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(data))
module_trait_padj <- matrix(
  p.adjust(as.vector(module_trait_pval), method = "BH"),
  nrow = nrow(module_trait_pval),
  dimnames = dimnames(module_trait_pval)
)

full_results <- left_join(
  as.data.frame(module_trait_cor) %>%
    tibble::rownames_to_column("module") %>%
    tidyr::pivot_longer(-module, names_to = "trait", values_to = "correlation"),
  as.data.frame(module_trait_padj) %>%
    tibble::rownames_to_column("module") %>%
    tidyr::pivot_longer(-module, names_to = "trait", values_to = "padj"),
  by = c("module", "trait")
)

write.csv(full_results,
          file.path(params$out_dir, "module_trait_correlations.csv"),
          row.names = FALSE)

# Heatmap
text_matrix <- paste0(signif(module_trait_cor, 2), "\n(",
                      signif(module_trait_padj, 1), ")")
dim(text_matrix) <- dim(module_trait_cor)

pdf(file.path(params$out_dir, "03_module_trait_heatmap.pdf"), width = 12, height = 8)
par(mar = c(6, 9, 3, 3))
labeledHeatmap(
  Matrix        = module_trait_cor,
  xLabels       = colnames(conds),
  yLabels       = rownames(module_trait_cor),
  ySymbols      = rownames(module_trait_cor),
  colorLabels   = FALSE,
  colors        = blueWhiteRed(50),
  textMatrix    = text_matrix,
  setStdMargins = FALSE,
  cex.text      = 0.5,
  zlim          = c(-1, 1),
  main          = "Module-Trait Relationships (r, adj. p-value)"
)
dev.off()


# Identify Significant Modules --------------------------------------------


sig_modules <- full_results %>%
  filter(abs(correlation) >= params$cor_threshold,
         padj             <  params$p_threshold) %>%
  mutate(module_color = sub("^ME", "", module)) %>%
  pull(module_color) %>%
  unique()

cat(sprintf("\nSignificant modules (n = %d): %s\n",
            length(sig_modules), paste(sig_modules, collapse = ", ")))

if (length(sig_modules) == 0)
  stop("No modules passed thresholds. Relax cor_threshold or p_threshold.")


# Build Combined Edge List & Hub Gene Table -------------------------------


probes    <- colnames(data)
all_edges <- list()
all_hubs  <- list()

for (mod in sig_modules) {
  
  cat(sprintf("\nProcessing module: %s\n", mod))
  
  in_module  <- module_colors == mod
  mod_probes <- probes[in_module]
  
  # Subset TOM to this module
  mod_TOM             <- (1 - dissTOM)[in_module, in_module]
  dimnames(mod_TOM)   <- list(mod_probes, mod_probes)
  
  # Extract upper triangle above threshold (undirected, no duplicates)
  upper_idx <- which(upper.tri(mod_TOM) & mod_TOM >= params$tom_threshold,
                     arr.ind = TRUE)
  
  if (nrow(upper_idx) == 0) {
    warning(sprintf("Module %s: no edges above threshold %.2f — skipping.",
                    mod, params$tom_threshold))
    next
  }
  
  # --- Edge list: source_gene | target_gene | correlation | module
  mod_edges <- data.frame(
    source_gene = mod_probes[upper_idx[, 1]],
    target_gene = mod_probes[upper_idx[, 2]],
    correlation = mod_TOM[upper_idx],
    module      = mod,
    stringsAsFactors = FALSE
  )
  all_edges[[mod]] <- mod_edges
  
  # --- Hub genes via igraph centrality metrics
  g <- graph_from_data_frame(
    d        = mod_edges %>% select(from = source_gene,
                                    to   = target_gene,
                                    weight = correlation),
    directed = FALSE,
    vertices = mod_probes
  )
  
  node_metrics <- data.frame(
    gene        = V(g)$name,
    module      = mod,
    degree      = degree(g),
    betweenness = betweenness(g, normalized = TRUE),
    closeness   = closeness(g,  normalized = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Hub = top hub_quantile by degree within this module
  hub_cutoff       <- quantile(node_metrics$degree, params$hub_quantile)
  all_hubs[[mod]]  <- node_metrics %>%
    filter(degree >= hub_cutoff) %>%
    arrange(desc(betweenness))
  
  cat(sprintf("  Genes: %d | Edges: %d | Hubs: %d\n",
              length(mod_probes), nrow(mod_edges), nrow(all_hubs[[mod]])))
}


# Save Outputs ------------------------------------------------------------


# Combined edge list
edge_list <- bind_rows(all_edges)
write.csv(edge_list,
          file.path(params$out_dir, "combined_edge_list.csv"),
          row.names = FALSE)

# Hub gene table
hub_table <- bind_rows(all_hubs)
write.csv(hub_table,
          file.path(params$out_dir, "hub_genes.csv"),
          row.names = FALSE)

# Export tan module gene list
tan_genes <- data.frame(
  gene   = probes[module_colors == "tan"],
  module = "tan"
)

write.csv(tan_genes,
          file.path(params$out_dir, "tan_module_genes.csv"),
          row.names = FALSE)

cat(sprintf("Tan module: %d genes exported\n", nrow(tan_genes)))


# Summary -----------------------------------------------------------------


cat("\n--- Output Summary ---\n")
cat(sprintf("combined_edge_list.csv : %d edges across %d modules\n",
            nrow(edge_list), length(all_edges)))
cat(sprintf("hub_genes.csv          : %d hub genes across %d modules\n",
            nrow(hub_table), length(all_hubs)))
cat(sprintf("Output directory       : %s\n", params$out_dir))

cat("\nEdges per module:\n")
edge_list %>% count(module, name = "n_edges") %>% print()

cat("\nHub genes per module:\n")
hub_table %>% count(module, name = "n_hubs") %>% print()