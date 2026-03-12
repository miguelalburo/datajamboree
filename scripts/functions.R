# 
# Functions to streamline the code in the main script
#


# Loading Libraries & Packages --------------------------------------------

# Clear Environment
rm(list = ls())


# Load libraries
library(tidyverse)


# Load data
metadata <- read.csv('../data/sample_sheet.csv') # sample metadata


# Processing Functions ----------------------------------------------------


# 'Omics Loader
# Convert gene counts / DIMS peak intensity to a row = sample, column = feature 
load_feature_matrix <- function(path, sep = ",") {

  df <- read.delim(path, sep = sep, row.names = 1) %>%
    t() %>% as.data.frame()
  
  return(df)
}


# Condition Filter
# Select Samples based on REF, Treatment, Site
condition_filter <- function(site = NULL, ref = NULL, desc = NULL) {
  
  samples <- metadata %>%
    filter(
      if (!is.null(site)) Site %in% site        else TRUE,
      if (!is.null(ref))  REF %in% ref          else TRUE,
      if (!is.null(desc)) Description %in% desc else TRUE
    ) %>%
    pull(SampleID)
  
  return(samples)
}


# Plotting Functions ------------------------------------------------------


# Factor plot for RGCCA
factor.plot <- function(y1, y2, conds){
  
  df <- data.frame(datablock1 = y1, datablock2 = y2, color = conds)
  
  p <- ggplot(df, aes(datablock1, datablock2, color = color)) + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(aes(alpha = 0.6, size = 1)) +          # Moved outside aes()
    guides(size = FALSE, alpha = FALSE) + # Remove size and alpha legends
    labs(subtitle = "Component 1 factor plot") +
    theme_bw() +
    theme(legend.position = "right")
  return(p)
}

