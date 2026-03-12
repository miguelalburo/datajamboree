

```{r}
# Load required packages
library(ggplot2)
library(factoextra)
library(pheatmap)
library(vegan)

# Load the helper function for metadata-data alignment
source("fileio.R")



## Positive mode

# Load positive-mode metabolomics data with matched metadata
pos.ds <- load.dataset(
  meta.file = "sample_sheet.csv",
  meta.sep = ",",
  data.file = "polar_pos_pqn_imputed_glog.csv",
  data.sep = ","
)

# Separate matrix and labels
pos.features <- pos.ds$data.matrix
pos.labels <- pos.ds$meta.data

# Make grouping variables explicit factors
pos.labels$REF <- factor(pos.labels$REF, levels = c("Control", "1x", "10x"))
pos.labels$Site <- factor(pos.labels$Site)

# Check the structure
str(pos.features)
str(pos.labels)

# Check whether any NA values remain after imputation
sum(is.na(pos.features))

# Identify constant columns
pos.const_cols <- apply(pos.features, 2, function(x) sd(x, na.rm = TRUE) == 0)

# Show constant feature names if present
colnames(pos.features)[pos.const_cols]

# Keep only columns with variation
pos.features.pca <- pos.features[, !pos.const_cols, drop = FALSE]

# Run PCA
pos.pca <- prcomp(pos.features.pca, center = TRUE, scale. = TRUE)

# Scree plot
fviz_eig(pos.pca)

# PCA by concentration group
fviz_pca_biplot(
  pos.pca,
  geom = "point",
  habillage = pos.labels$REF,
  addEllipses = TRUE,
  invisible = "var"
) +
  ggtitle("Positive mode PCA (coloured by REF)")

# PCA by site
fviz_pca_biplot(
  pos.pca,
  geom = "point",
  habillage = pos.labels$Site,
  addEllipses = TRUE,
  invisible = "var"
) +
  ggtitle("Positive mode PCA (coloured by Site)")



## Positive mode sample distance heatmap

# Scale features before distance calculation
pos.scaled <- scale(pos.features.pca)

# Calculate Euclidean distance between samples
pos.dist <- dist(pos.scaled)

# Prepare sample annotations for the heatmap
pos.annotation <- pos.labels[, c("REF", "Site"), drop = FALSE]

# Draw the sample distance heatmap
pheatmap(
  as.matrix(pos.dist),
  annotation_row = pos.annotation,
  annotation_col = pos.annotation,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Positive mode sample distance heatmap"
)



## Positive mode PERMANOVA

# Keep treated samples only for factorial testing
pos.labels.tr <- subset(pos.labels, REF != "Control")
pos.features.tr <- pos.features.pca[rownames(pos.labels.tr), , drop = FALSE]

# Scale treated data before distance calculation
pos.scaled.tr <- scale(pos.features.tr)

# Run PERMANOVA on treated samples only
# This tests Site, REF (1x vs 10x), and Site:REF interaction
# Control is excluded here because it is not site-matched
pos.adonis <- adonis2(
  dist(pos.scaled.tr) ~ Site * REF,
  data = pos.labels.tr
)

# Show PERMANOVA results
pos.adonis



## Negative mode

# Load negative-mode metabolomics data with matched metadata
neg.ds <- load.dataset(
  meta.file = "sample_sheet.csv",
  meta.sep = ",",
  data.file = "polar_neg_pqn_imputed_glog.csv",
  data.sep = ","
)

# Separate matrix and labels
neg.features <- neg.ds$data.matrix
neg.labels <- neg.ds$meta.data

# Make grouping variables explicit factors
neg.labels$REF <- factor(neg.labels$REF, levels = c("Control", "1x", "10x"))
neg.labels$Site <- factor(neg.labels$Site)

# Check the structure
str(neg.features)
str(neg.labels)

# Check whether any NA values remain after imputation
sum(is.na(neg.features))

# Identify constant columns
neg.const_cols <- apply(neg.features, 2, function(x) sd(x, na.rm = TRUE) == 0)

# Show constant feature names if present
colnames(neg.features)[neg.const_cols]

# Keep only columns with variation
neg.features.pca <- neg.features[, !neg.const_cols, drop = FALSE]

# Run PCA
neg.pca <- prcomp(neg.features.pca, center = TRUE, scale. = TRUE)

# Scree plot
fviz_eig(neg.pca)

# PCA by concentration group
fviz_pca_biplot(
  neg.pca,
  geom = "point",
  habillage = neg.labels$REF,
  addEllipses = TRUE,
  invisible = "var"
) +
  ggtitle("Negative mode PCA (coloured by REF)")

# PCA by site
fviz_pca_biplot(
  neg.pca,
  geom = "point",
  habillage = neg.labels$Site,
  addEllipses = TRUE,
  invisible = "var"
) +
  ggtitle("Negative mode PCA (coloured by Site)")



## Negative mode sample distance heatmap

# Scale features before distance calculation
neg.scaled <- scale(neg.features.pca)

# Calculate Euclidean distance between samples
neg.dist <- dist(neg.scaled)

# Prepare sample annotations for the heatmap
neg.annotation <- neg.labels[, c("REF", "Site"), drop = FALSE]

# Draw the sample distance heatmap
pheatmap(
  as.matrix(neg.dist),
  annotation_row = neg.annotation,
  annotation_col = neg.annotation,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Negative mode sample distance heatmap"
)



## Negative mode PERMANOVA

# Keep treated samples only for factorial testing
neg.labels.tr <- subset(neg.labels, REF != "Control")
neg.features.tr <- neg.features.pca[rownames(neg.labels.tr), , drop = FALSE]

# Scale treated data before distance calculation
neg.scaled.tr <- scale(neg.features.tr)

# Run PERMANOVA on treated samples only
# This tests Site, REF (1x vs 10x), and Site:REF interaction
# Control is excluded here because it is not site-matched
neg.adonis <- adonis2(
  dist(neg.scaled.tr) ~ Site * REF,
  data = neg.labels.tr
)

# Show PERMANOVA results
neg.adonis
```