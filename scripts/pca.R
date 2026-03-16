#
# Principal Component Analysis (PCA)
#


# Setup -------------------------------------------------------------------

# Essential
source("functions.R")

library(dplyr)
library(ggplot2)
library(ggrepel)    # non-overlapping labels: install.packages("ggrepel")


# Functions ---------------------------------------------------------------


#' Run PCA on a samples x genes expression matrix
#'
#' @param data      Numeric matrix (samples x features)
#' @param meta      Data frame of sample metadata (rows must match data rows)
#' @param scale     Logical; scale features to unit variance (default TRUE)
#' @param center    Logical; center features to zero mean (default TRUE)
#'
#' @return A list:
#'   $scores     — PC scores for each sample (data frame, includes metadata)
#'   $variance   — Variance explained per PC (data frame)
#'   $loadings   — Feature loadings (data frame)
#'   $prcomp     — Raw prcomp object for further use

run_pca <- function(data, meta, scale = TRUE, center = TRUE) {
  
  # Sanity checks
  if (!identical(rownames(data), rownames(meta)))
    stop("Row names of data and meta must match.")
  if (any(is.na(data)))
    stop("data contains NA values. Impute or remove before running PCA.")
  if (any(apply(data, 2, var) == 0))
    message("Removing zero-variance features before PCA.")
  
  # Remove zero-variance features (PCA will fail on them)
  data <- data[, apply(data, 2, var) != 0]
  
  pca <- prcomp(data, center = center, scale. = scale)
  
  # Variance explained
  var_explained <- data.frame(
    PC              = paste0("PC", seq_along(pca$sdev)),
    variance        = pca$sdev^2,
    pct_variance    = (pca$sdev^2 / sum(pca$sdev^2)) * 100
  ) %>%
    mutate(cumulative_pct = cumsum(pct_variance))
  
  # Sample scores joined with metadata
  scores <- as.data.frame(pca$x) %>%
    tibble::rownames_to_column("sample_id") %>%
    bind_cols(meta)
  
  # Feature loadings
  loadings <- as.data.frame(pca$rotation) %>%
    tibble::rownames_to_column("feature")
  
  list(
    scores   = scores,
    variance = var_explained,
    loadings = loadings,
    prcomp   = pca
  )
}


#' Scree plot — variance explained per PC
#'
#' @param variance    $variance data frame from run_pca()
#' @param n_pcs       Number of PCs to show (default 20)

plot_scree <- function(variance, n_pcs = 20) {
  
  plot_data <- variance %>%
    slice_head(n = n_pcs) %>%
    mutate(PC = factor(PC, levels = PC))
  
  ggplot(plot_data, aes(x = PC, y = pct_variance)) +
    geom_col(fill = "#457B9D", alpha = 0.85) +
    geom_line(aes(group = 1), colour = "#E63946", linewidth = 0.7) +
    geom_point(colour = "#E63946", size = 2) +
    geom_text(aes(label = sprintf("%.1f%%", pct_variance)),
              vjust = -0.5, size = 3, colour = "grey30") +
    labs(
      title = "Scree Plot",
      x     = NULL,
      y     = "Variance Explained (%)"
    ) +
    theme_classic() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      plot.title   = element_text(face = "bold")
    )
}


#' PCA scatter plot — samples coloured by a metadata variable
#'
#' @param scores      $scores data frame from run_pca()
#' @param variance    $variance data frame from run_pca()
#' @param colour_by   Column name in scores to colour points by
#' @param shape_by    Column name in scores to shape points by (optional)
#' @param label_by    Column name in scores to label points by (optional)
#' @param pc_x        PC to plot on x-axis (default 1)
#' @param pc_y        PC to plot on y-axis (default 2)
#' @param point_size  Point size (default 3)
#' @param palette     Named or unnamed vector of colours (optional)

plot_pca <- function(scores,
                     variance,
                     colour_by,
                     title = NULL,
                     shape_by   = NULL,
                     label_by   = NULL,
                     pc_x       = 1,
                     pc_y       = 2,
                     point_size = 3,
                     palette    = NULL) {
  
  x_col  <- paste0("PC", pc_x)
  y_col  <- paste0("PC", pc_y)
  x_pct  <- round(variance$pct_variance[pc_x], 1)
  y_pct  <- round(variance$pct_variance[pc_y], 1)
  
  aes_mapping <- aes(
    x      = .data[[x_col]],
    y      = .data[[y_col]],
    colour = .data[[colour_by]]
  )
  
  p <- ggplot(scores, aes_mapping) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey80") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey80") +
    geom_point(
      aes(shape = if (!is.null(shape_by)) .data[[shape_by]] else NULL),
      size  = point_size,
      alpha = 0.85
    ) +
    labs(
      title = title,
      x      = sprintf("%s  (%.1f%%)", x_col, x_pct),
      y      = sprintf("%s  (%.1f%%)", y_col, y_pct),
      colour = colour_by,
      shape  = shape_by
    ) +
    theme_classic() +
    theme(
      legend.position = "right",
      plot.title      = element_text(face = "bold")
    )
  
  # Optional sample labels (non-overlapping via ggrepel)
  if (!is.null(label_by)) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = .data[[label_by]]),
      size          = 3,
      max.overlaps  = 20,
      colour        = "grey20",
      segment.colour = "grey70"
    )
  }
  
  # Custom palette
  if (!is.null(palette)) {
    p <- p + scale_colour_manual(values = palette)
  } else {
    p <- p + scale_colour_brewer(palette = "Paired")
  }
  
  p
}


#' Top N feature loadings for a given PC — useful for interpreting drivers
#'
#' @param loadings    $loadings data frame from run_pca()
#' @param pc          PC number to inspect (default 1)
#' @param n           Number of top features to return (default 20)

plot_loadings <- function(loadings, pc = 1, n = 20) {
  
  pc_col <- paste0("PC", pc)
  
  plot_data <- loadings %>%
    select(feature, loading = all_of(pc_col)) %>%
    slice_max(abs(loading), n = n) %>%
    arrange(loading) %>%
    mutate(feature   = factor(feature, levels = feature),
           direction = if_else(loading >= 0, "Positive", "Negative"))
  
  ggplot(plot_data, aes(x = loading, y = feature, fill = direction)) +
    geom_col(alpha = 0.85) +
    geom_vline(xintercept = 0, colour = "grey40") +
    scale_fill_manual(values = c("Positive" = "#2DC653", "Negative" = "#E63946")) +
    labs(
      title = sprintf("Top %d Loadings — PC%d", n, pc),
      x     = "Loading",
      y     = NULL,
      fill  = NULL
    ) +
    theme_classic() +
    theme(
      legend.position = "top",
      plot.title      = element_text(face = "bold")
    )
}


# Datasets ----------------------------------------------------------------


# RNAseq Gene Counts
rna.ds <- load.dataset(
  meta.file = '../data/sample_sheet.csv', meta.sep = ',',
  data.file = '../data/rna_vst_counts.csv', data.sep = ','
)

# DIMS Peak Intensity
dims.ds <- load.dataset(
  meta.file = '../data/sample_sheet.csv', meta.sep = ',',
  data.file = '../data/polar_pos_pqn_imputed_glog.csv', data.sep = ','
)

# Combined List
datasets <- list(transcriptomics = rna.ds, metabolomics = dims.ds)


# Run PCA -----------------------------------------------------------------


pca_results <- lapply(datasets, function(ds) {
  
  data <- ds[["data.matrix"]]
  
  meta <- ds[["meta.data"]]
  
  run_pca(data, meta, scale = TRUE, center = TRUE)
})


# Plots & Export ----------------------------------------------------------


out_dir <- "../results/pca"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (name in names(datasets)) {
  
  pca   <- pca_results[[name]]
  label <- datasets[[name]]$label
  
  cat(sprintf("\n--- %s ---\n", label))
  print(pca$variance %>% filter(PC %in% paste0("PC", 1:10)))
  
  # Scree plot
  p_scree <- plot_scree(pca$variance, n_pcs = 15) +
    ggtitle(sprintf("Scree Plot — %s", label))
  
  # PC1 vs PC2
  p_pca12 <- plot_pca(
    scores    = pca$scores,
    variance  = pca$variance,
    colour_by = "Site",     # change to your metadata column
    shape_by  = "REF",         # set NULL if not needed
    label_by  = NULL,     # set NULL to hide labels
    title = str_to_sentence(name),
    pc_x      = 1,
    pc_y      = 2
  ) + ggtitle(sprintf("PCA — %s", label))
  
  # Top 20 loadings for PC1
  p_load <- plot_loadings(pca$loadings, pc = 1, n = 20) +
    ggtitle(sprintf("Top 20 Loadings PC1 — %s", label))
  
  # Save plots
  ggsave(file.path(out_dir, sprintf("%s_scree.png",       name)), p_scree,  width = 8, height = 5, dpi = 300)
  ggsave(file.path(out_dir, sprintf("%s_pca_pc1_pc2.png", name)), p_pca12,  width = 8, height = 6, dpi = 300)
  ggsave(file.path(out_dir, sprintf("%s_loadings_pc1.png",name)), p_load,   width = 7, height = 8, dpi = 300)
  
  # Save tables
  write.csv(pca$scores,   file.path(out_dir, sprintf("%s_scores.csv",   name)), row.names = FALSE)
  write.csv(pca$variance, file.path(out_dir, sprintf("%s_variance.csv", name)), row.names = FALSE)
  write.csv(pca$loadings, file.path(out_dir, sprintf("%s_loadings.csv", name)), row.names = FALSE)
  
  cat(sprintf("%s outputs written to: %s\n", label, out_dir))
}

cat("\nAll datasets complete.\n")