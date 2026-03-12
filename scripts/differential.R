# Differential expression analysis using DESeq2
rm(list=ls())

# Load necessary libraries
library(DESeq2)

# Load data
counts <- read.csv('../data/rna_vst_counts.csv', row.names = 1)
metadata <- read.csv('../data/sample_sheet.csv', row.names = 1)

# Ensure that the sample names in colData match the column names in countData
metadata <- metadata[match(colnames(counts), rownames(metadata)), ]

# Formatting
counts <- round(counts)
metadata$Site <- factor(metadata$Site)
metadata$REF <- factor(metadata$REF)

# Build DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = counts,   # DESeq2 needs integers
  colData   = metadata,
  design    = ~ Site
)

# Run LRT — drop 'site' in the reduced model
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
res_lrt  <- results(dds_lrt)

# How many genes respond across sites?
sum(res_lrt$padj < 0.05, na.rm = TRUE)
