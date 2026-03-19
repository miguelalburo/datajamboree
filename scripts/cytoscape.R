#
# Script for visualising Gene Co-expression Networks in Cytoscape
# 


# Setup -------------------------------------------------------------------

source("functions.R")

library(RCy3)
library(RColorBrewer)

NETWORK_DIR <- "../results/wgcna/"
EXPORT_DIR  <- "../results/wgcna/figures/"

dir.create(EXPORT_DIR, showWarnings = FALSE, recursive = TRUE)

cytoscapePing()


# Load combined edge list -------------------------------------------------

# Expected columns: source | target | weight | module
edge_list <- read.csv(file.path(NETWORK_DIR, "combined_edge_list.csv"),
                      stringsAsFactors = FALSE) %>% 
  rename(source = 1, target = 2, weight = 3) %>%
  filter(weight > 0.35)


stopifnot(all(c("source", "target", "weight", "module") %in% names(edge_list)))

modules <- unique(edge_list$module)
message(sprintf("Found %d module(s): %s", length(modules), paste(modules, collapse = ", ")))


# Helper: resolve module colour to hex ------------------------------------

module_to_hex <- function(m) {
  # If the module name is itself a recognised R colour (common in WGCNA),
  # convert it directly; otherwise fall back to a Set2 palette colour.
  if (m %in% grDevices::colours()) {
    grDevices::rgb(t(grDevices::col2rgb(m)), maxColorValue = 255)
  } else {
    pal <- RColorBrewer::brewer.pal(max(3, length(modules)), "Set2")
    pal[which(modules == m) %% length(pal) + 1]
  }
}


# Batch loop over modules -------------------------------------------------

for (mod in modules) {
  
  message("\nProcessing module: ", mod)
  
  # --- Subset edges for this module --------------------------------------
  
  edges <- edge_list[edge_list$module == mod, c("source", "target", "weight")]
  
  if (nrow(edges) == 0) {
    message("  Skipping — no edges found.")
    next
  }
  
  
  # --- Build node table from edge endpoints ------------------------------
  
  all_ids <- unique(c(edges$source, edges$target))
  
  degree_tbl <- dplyr::bind_rows(
    dplyr::select(edges, id = source),
    dplyr::select(edges, id = target)
  ) %>%
    dplyr::count(id, name = "degree")
  
  connectivity_tbl <- dplyr::bind_rows(
    dplyr::select(edges, id = source, weight),
    dplyr::select(edges, id = target, weight)
  ) %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(connectivity = sum(weight, na.rm = TRUE), .groups = "drop")
  
  nodes <- data.frame(id = all_ids, module = mod, stringsAsFactors = FALSE) %>%
    dplyr::left_join(degree_tbl,       by = "id") %>%
    dplyr::left_join(connectivity_tbl, by = "id") %>%
    dplyr::mutate(
      degree       = dplyr::coalesce(degree, 0L),
      connectivity = dplyr::coalesce(connectivity, 0),
      is_hub       = connectivity >= quantile(connectivity, 0.95)
    )
  
  
  # --- Resolve colour for this module ------------------------------------
  
  hex_colour <- module_to_hex(mod)
  
  
  # --- Create network in Cytoscape ---------------------------------------
  
  title <- paste0("Module_", mod)
  
  net_suid <- createNetworkFromDataFrames(
    nodes      = nodes,
    edges      = edges,
    title      = title,
    collection = "WGCNA"
  )
  
  
  # --- Visual style ------------------------------------------------------
  
  style_name <- paste0("CoExp_", title)
  createVisualStyle(style_name)
  setVisualStyle(style_name)
  
  setBackgroundColorDefault("#1A1A2E", style.name = style_name)
  setNodeLabelDefault("",             style.name = style_name)
  
  # All nodes in this module share the same module colour
  setNodeColorDefault(hex_colour, style.name = style_name)
  
  # Node size by connectivity
  setNodeSizeMapping(
    table.column        = "connectivity",
    table.column.values = c(min(nodes$connectivity), max(nodes$connectivity)),
    sizes               = c(15, 70),
    mapping.type        = "c",
    style.name          = style_name
  )
  
  # Hub nodes: white border + diamond shape
  setNodeBorderColorMapping(
    table.column        = "is_hub",
    table.column.values = c("TRUE", "FALSE"),
    colors              = c("#FFFFFF", "#333333"),
    mapping.type        = "d",
    style.name          = style_name
  )
  setNodeBorderWidthMapping(
    table.column        = "is_hub",
    table.column.values = c("TRUE", "FALSE"),
    widths              = c(3, 0.5),
    mapping.type        = "d",
    style.name          = style_name
  )
  setNodeShapeMapping(
    table.column        = "is_hub",
    table.column.values = c("TRUE", "FALSE"),
    shapes              = c("DIAMOND", "ELLIPSE"),
    style.name          = style_name
  )
  
  # Edge width and opacity by weight
  setEdgeLineWidthMapping(
    table.column        = "weight",
    table.column.values = c(min(edges$weight), max(edges$weight)),
    widths              = c(0.5, 6),
    mapping.type        = "c",
    style.name          = style_name
  )
  setEdgeOpacityMapping(
    table.column        = "weight",
    table.column.values = c(min(edges$weight), max(edges$weight)),
    opacities           = c(60, 220),
    mapping.type        = "c",
    style.name          = style_name
  )
  setEdgeColorDefault("#7A7A8C",   style.name = style_name)
  setEdgeTargetArrowShapeDefault("NONE", style.name = style_name)
  setEdgeSourceArrowShapeDefault("NONE", style.name = style_name)
  
  
  # --- Layout & export ---------------------------------------------------
  
  layoutNetwork(
    "force-directed defaultSpringCoefficient=0.00004 defaultSpringLength=60",
    network = net_suid
  )
  fitContent()
  
  exportImage(
    filename      = file.path(EXPORT_DIR, paste0(title, ".png")),
    type          = "PNG",
    resolution    = 300,
    overwriteFile = TRUE
  )
  exportNetwork(
    filename      = file.path(EXPORT_DIR, paste0(title, ".graphml")),
    type          = "graphml",
    overwriteFile = TRUE
  )
  
  message("  Done — ", nrow(nodes), " nodes, ", nrow(edges), " edges")
}

message("\nAll modules processed.")