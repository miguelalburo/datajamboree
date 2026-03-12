# Preamble ---------------------------------------------------------------------

# Formulate our genes-of-interest
# Approach with DEGs.
# Set working directory to source location of this script

rm(list=ls())

# Load necessary libraries
library(tidyverse)
library(caret)
library(limma)
library(pheatmap)

# Load data
counts_t <- read.csv('../data/rna_vst_counts.csv', row.names = 1)
metadata <- read.csv('../data/sample_sheet.csv')

# Pre-processing ----------------------------------------------------------

# Classic expression matrix format
counts <- counts_t %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "sampleid")


# PCA ---------------------------------------------------------------------

# PCA
pca_res <- counts %>%
  select(-sampleid) %>%
  prcomp(center = TRUE, scale. = TRUE)

summary(pca_res)

# Scree plot
variance <- data.frame(
  pc = c(1:length(pca_res$sdev)),
  var = pca_res$sdev^2 / sum(pca_res$sdev^2)
) %>% head(10)

p.scree <- ggplot(variance, aes(x = pc, y = var)) +
  geom_line() +
  geom_point(color = "red", show.legend = F) +
  labs(title = "Scree Plot",
       x = "Principal Component",
       y = "Proportion of Variance Explained") +
  scale_x_continuous(breaks = 1:10) +
  theme_bw()

print(p.scree)

# PC1-PC2 Biplot
pca_df <- data.frame(
  sampleid = counts$sampleid,
  site = metadata$Site[match(counts$sampleid, metadata$SampleID)],
  concentration = metadata$REF[match(counts$sampleid, metadata$SampleID)],
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2]
)

p.biplot <- ggplot(pca_df, 
                   aes(x = PC1, y = PC2,
                       color = site, shape = concentration)) +
  geom_point() +
  labs(title = "PCA Biplot",
       x = paste("PC1",round(variance$var[1]*100, 2), "%"),
       y = paste("PC2",round(variance$var[2]*100, 2), "%")) +
  theme_bw()

print(p.biplot)

# Differential Expression ------------------------------------------------------

# Model design for limma
design <- model.matrix(~ Site, data = metadata)

# Fit linear model and compute statistics
fit <- lmFit(counts_t, design)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Global ANOVA-style test
anova_res <- topTable(fit,
                      coef = 2:ncol(design),
                      number = Inf,
                      sort.by = "F")


# max log2 fold change across conditions
max_fc <- apply(counts_t, 1, function(x) max(x) - min(x))

anova_res$max_log2FC <- max_fc[rownames(anova_res)]


# DEGs ranked by adjusted p-value and max log2 fold change
genes_ranked <- anova_res %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(adj.P.Val) %>%
  rownames_to_column(var = "geneid")

top50 <- genes_ranked$geneid[1:50] # top 50 DEGs


# Top DEG expression across sites 
df <- counts %>%
  mutate(site = match(sampleid, metadata$SampleID) %>% 
           metadata$Site[.]) %>%
  select(sampleid, site,
         all_of(genes_ranked$geneid[1:10])) %>%
  pivot_longer(cols = -c(sampleid, site),
               names_to = "gene", values_to = "expression")

p.boxplot <- df %>%
  ggplot(aes(x = site, y = expression, fill = site)) +
  geom_boxplot() +
  facet_wrap(~ gene, scales = "free_y") +
  labs(title = "Top 10 DEGs Expression Across Sites",
       x = "Site",
       y = "Expression (VST)") +
  theme_bw() +
  theme(legend.position = "none")
print(p.boxplot)

# Heatmap of top DEGs ------------------------------------------------------

counts %>% 
  column_to_rownames(var = "sampleid") %>%
  select(all_of(top50)) %>% 
  pheatmap(show_rownames = F, show_colnames = F)
  

# Outputs -----------------------------------------------------------------

# DEGs
write.csv(genes_ranked,"../results/degs.csv")

