#
# Script for visualising Gene Co-expression Networks in Cytoscape
# 
#


# Setup -------------------------------------------------------------------

# Essentials
source("functions.R")

# Specifics
library(RCy3)
library(RColorBrewer)


# Configs
NETWORK_DIR  <- "../results/wgcna/"
EXPORT_DIR   <- "../results/wgcna/figures/"

dir.create(EXPORT_DIR, showWarnings = FALSE, recursive = TRUE)

# Ping Cytoscape
cytoscapePing()


# Discover file pairs -----------------------------------------------------


node_files <- list.files(NETWORK_DIR, pattern = "_nodes\\.txt$", full.names = TRUE)
edge_files <- sub("_nodes\\.txt$", "_edges.txt", node_files)

# Drop pairs where the edge file is missing
keep       <- file.exists(edge_files)
node_files <- node_files[keep]
edge_files <- edge_files[keep]

message(sprintf("Found %d network pair(s)", length(node_files)))


# Batch loop --------------------------------------------------------------


for (i in seq_along(node_files)) {
  
  nf    <- node_files[i]
  ef    <- edge_files[i]
  title <- sub("_nodes\\.txt$", "", basename(nf))
  
  message("\nProcessing: ", title)
  
  
  # Load data -------------------------------------------------------------
  
  nodes <- read.delim(nf, header = TRUE, sep = "\t") %>%
    as.data.frame() %>%
    select(-any_of("altName")) %>%
    rename(id = "nodeName", group = "nodeAttr.nodesPresent...")
  
  edges <- read.delim(ef, header = TRUE, sep = "\t") %>%
    as.data.frame() %>%
    select(-any_of(c("fromAltName", "toAltName"))) %>%
    rename(source = "fromNode", target = "toNode")
  
  
  # Node metrics ----------------------------------------------------------
  
  degree_tbl <- bind_rows(
    edges %>% select(id = source),
    edges %>% select(id = target)
  ) %>%
    count(id, name = "degree")
  
  connectivity_tbl <- bind_rows(
    edges %>% select(id = source, weight),
    edges %>% select(id = target, weight)
  ) %>%
    group_by(id) %>%
    summarise(connectivity = sum(weight, na.rm = TRUE), .groups = "drop")
  
  nodes <- nodes %>%
    left_join(degree_tbl,       by = "id") %>%
    left_join(connectivity_tbl, by = "id") %>%
    mutate(
      degree       = replace(degree,       is.na(degree),       0),
      connectivity = replace(connectivity, is.na(connectivity), 0),
      is_hub       = connectivity >= quantile(connectivity, 0.95)
    )
  
  
  # Module colours --------------------------------------------------------
  
  modules <- unique(nodes$group)
  
  module_colors <- setNames(
    sapply(modules, function(m) {
      if (m %in% grDevices::colours()) {
        grDevices::rgb(t(grDevices::col2rgb(m)), maxColorValue = 255)
      } else {
        brewer.pal(max(3, length(modules)), "Set2")[which(modules == m)]
      }
    }),
    modules
  )
  
  
  # Create network --------------------------------------------------------
  
  net_suid <- createNetworkFromDataFrames(
    nodes      = nodes,
    edges      = edges,
    title      = title,
    collection = "WGCNA"
  )
  
  
  # Visual style ----------------------------------------------------------
  
  style_name <- paste0("CoExp_", title)
  createVisualStyle(style_name)
  setVisualStyle(style_name)
  
  setBackgroundColorDefault("#1A1A2E", style.name = style_name)
  
  # No labels on any nodes
  setNodeLabelDefault("", style.name = style_name)
  
  # Node colour by module
  setNodeColorMapping(
    table.column        = "group",
    table.column.values = modules,
    colors              = unname(module_colors),
    mapping.type        = "d",
    style.name          = style_name
  )
  
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
  setEdgeColorDefault("#7A7A8C", style.name = style_name)
  setEdgeTargetArrowShapeDefault("NONE", style.name = style_name)
  setEdgeSourceArrowShapeDefault("NONE", style.name = style_name)
  
  
  # Layout ----------------------------------------------------------------
  
  layoutNetwork(
    "force-directed defaultSpringCoefficient=0.00004 defaultSpringLength=60",
    network = net_suid
  )
  fitContent()
  
  
  # Export ----------------------------------------------------------------
  
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

message("\nAll networks processed.")