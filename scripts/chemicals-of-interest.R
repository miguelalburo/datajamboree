# Preamble ---------------------------------------------------------------------
rm(list=ls())

# Load necessary libraries
library(tidyverse)

# Load data
chemicals_t <- read.delim('../data/water_chemicals.tsv') # chemical concentrations
counts_t <- read.csv('../data/rna_vst_counts.csv', row.names = 1) # gene counts
degs <- read.csv('../results/degs.csv') # differentially expressed genes
metadata <- read.csv('../data/sample_sheet.csv') # sample metadata

# Pre-processing ---------------------------------------------------------------

# Transposing counts to feature matrix
counts <- counts_t %>%
  t() %>% as.data.frame()


# Sample Filtering: Remove chemicals that are mostly zeros across sites
zeros_threshold = 0.25 # proportion of zeros

chemicals_t <- chemicals_t %>%
  filter(rowSums(across(where(is.numeric), ~ .x == 0)) <= zeros_threshold * 12)


# Add control site to chemicals
chemicals_t <- chemicals_t %>%
  mutate(Control = 0)


# # Long format
# long_chemicals <- chemicals %>%
#   pivot_longer(cols = where(is.numeric),
#                names_to = "Site",
#                values_to = "Concentration")


# Transposing chemicals to feature matrix
chemicals_short <- chemicals_t %>%
  column_to_rownames(var = "CAS") %>%
  select(-ChemName) %>% t() %>% as.data.frame() %>%
  rownames_to_column(var = "Site")

chemicals <-left_join(metadata,chemicals_short, by = "Site") %>%
  select(-REF, -Description, -Site) %>%
  column_to_rownames(var = "SampleID")


# Optional: add noise to avoid zero variance
chemicals_noise <- chemicals %>%
  mutate(across(where(is.numeric), ~ .x + rnorm(n(), mean = 0, sd = 1e-6)))


# REF filter
samples_to_keep <- metadata %>%
  filter(REF %in% c("1x", "Control")) %>%
  pull(SampleID)


# Pivot both to long & join
gene_long <- counts %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "gene_id", values_to = "counts")
chem_long <- chemicals %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "cas_id", values_to = "concentration")
combined <- gene_long %>%
  left_join(chem_long, by = "SampleID")


# Linear Gene ~ Chemical eQT-C  ------------------------------------------------

gene_mat <- counts[samples_to_keep, ] %>% as.matrix()
chem_mat <- chemicals_noise[samples_to_keep, ] %>% as.matrix()

fit_all_pairs <- function(Y, X) {
  n <- nrow(Y)
  
  results <- map_dfr(colnames(X), function(chem) {
    x     <- X[, chem]
    x_std <- cbind(1, x)                         # design matrix with intercept
    
    # Coefficients for all genes at once
    coefs <- solve(t(x_std) %*% x_std) %*% t(x_std) %*% Y
    
    # Residuals and standard errors for all genes
    fitted   <- x_std %*% coefs
    residuals <- Y - fitted
    rss      <- colSums(residuals^2)
    se       <- sqrt(rss / (n - 2) * solve(t(x_std) %*% x_std)[2, 2])
    
    # t-stats and p-values
    t_stat  <- coefs[2, ] / se
    p_value <- 2 * pt(abs(t_stat), df = n - 2, lower.tail = FALSE)
    
    tibble(
      cas_id      = chem,
      gene_id     = colnames(Y),
      estimate    = coefs[2, ],
      std.error   = se,
      statistic   = t_stat,
      p.value     = p_value
    )
  })
  
  results %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
    arrange(p.adj)
}

results <- fit_all_pairs(gene_mat, chem_mat)
tophits <- results %>% filter(p.adj < 0.05)


# Results ----------------------------------------------------------------------

# Numbers
tophits %>% nrow() # number of significant gene-chemical pairs
tophits$cas_id %>% unique() %>% length() # number of chemicals with significant associations
tophits$gene_id %>% unique() %>% length() # number of chemicals with significant associations


# Table
table_chemicals <- tophits$cas_id %>% table() %>% sort(decreasing=T)
table_genes <- tophits$gene_id %>% table() %>% sort(decreasing=T)



# Plots -------------------------------------------------------------------

# Scatterplot a gene vs a chemical
plot_gene_chem <- function(combined, gene_a, chem_b) {
  combined %>%
    filter(gene_id == gene_a, cas_id == chem_b) %>%
    ggplot(aes(x = concentration, y = counts)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, colour = "firebrick") +
    labs(title = paste(gene_a, "~", chem_b), x = "Concentration", y = "Counts")
}

geneA <- tophits$gene_id[2] # top gene
chemB <- tophits$cas_id[2] # top chemical
plot_gene_chem(combined, geneA, chemB)

# Number of gene-chemical associations by genes & chemicals
p.eqtc_gene <- table_genes %>%
  as.data.frame() %>%
  ggplot(aes(x = Freq)) +
  geom_histogram() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "Histogram of Significant Associations per Gene",
       x = "Number of Associated Chemicals per Gene",
       y = "Number of Genes")
p.eqtc_chem <- table_chemicals %>%
  as.data.frame() %>%
  ggplot(aes(x = Freq)) +
  geom_histogram() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "Histogram of Significant Associations per Chemical",
       x = "Number of Associated Genes per Chemical",
       y = "Number of Chemicals")