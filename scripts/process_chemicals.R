#
# Process Chemical Data
#
# - Filters Chemicals
# - Creates a sample x chemical matrix
# - Saves the matrix for later use in RGCCA
#

# Setup ---------------------------------------------------------------------


# Essentials
source("functions.R")


# Specific libraries
# library(tidyverse)

# Load data
chemicals_t <- read.delim('../data/water_chemicals.tsv') # chemicals
metadata <- read.csv('../data/sample_sheet.csv') # sample metadata

# Processing ---------------------------------------------------------------

# Zeros Filtering: Remove chemicals that are mostly zeros across sites
zeros_threshold = 0.25 # proportion of zeros

chemicals_t <- chemicals_t %>%
  filter(rowSums(across(where(is.numeric), ~ .x == 0)) <= zeros_threshold * 12)


# Add control site to chemicals
chemicals_t <- chemicals_t %>%
  mutate(Control = 0)


# Create CAS-ChemName mapping dataframe
cas_to_chemname <- chemicals_t %>%
  select(CAS, ChemName)


# Trasform to Classic Feature Matrix -------------------------------------------


# Transposing chemicals for rows = samples, columns = chemicals
chemicals_short <- chemicals_t %>%
  column_to_rownames(var = "CAS") %>%
  select(-ChemName) %>% t() %>% as.data.frame() %>%
  rownames_to_column(var = "Site")


# Site-to-SampleID Mapping (One-to-Many)
chemicals <-left_join(metadata,chemicals_short, by = "Site") %>%
  select(-REF, -Description, -Site) %>%
  column_to_rownames(var = "SampleID")


# Optional: add noise to avoid zero variance
chemicals_noise <- chemicals %>%
  mutate(across(where(is.numeric), ~ .x + rnorm(n(), mean = 0, sd = 1e-6)))


# Retranspose
chemicals <- chemicals %>%
  t() %>% as.data.frame()
chemicals_noise <- chemicals_noise %>%
  t() %>% as.data.frame()

# Save & Write ------------------------------------------------------------

# Imputated Concentrations Matrix
write.csv(chemicals, '../data/chem_conc.csv')
write.csv(chemicals_noise, '../data/chem_conc_noise.csv') # Noise


# CAS-ChemName Mapping
write.csv(cas_to_chemname, '../annotations/cas_to_chemname.csv')
