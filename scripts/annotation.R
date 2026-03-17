#
# Data and Function for Annotations
#



# Load Data ---------------------------------------------------------------

# Genes
hsa <- read.delim(file = "../annotations/rna_dma_to_hsa_mappings.tsv",
                  sep = "\t", header = F) %>%
  rename(daphnia = 1, ortholog = 2) %>% filter(ortholog != "")
  
dme <- read.delim(file = "../annotations/rna_dma_to_dme_mappings.tsv",
                  sep = "\t", header = F) %>% 
  rename(daphnia = 1, ortholog = 2) %>% filter(ortholog != "")

# Metabolites
pos.kegg <- read.delim(file = "../annotations/polar_pos_pkl_to_kegg_annotations.tsv",
                       sep = "\t", header = F) %>% 
  rename(mz = 1, metabolite = 3) %>% select(mz, metabolite) %>% 
  filter(metabolite != "")

neg.kegg <- read.delim(file = "../annotations/polar_neg_pkl_to_kegg_annotations.tsv",
                       sep = "\t", header = F) %>% 
  rename(mz = 1, metabolite = 3) %>% select(mz, metabolite) %>% 
  filter(metabolite != "")




# Mapping Functions -------------------------------------------------------

# Gene to Human Ensembl
daphnia_to_ortholog <- function(genes, mapping) {
  
  # Unique genes
  genes <- unique(genes)
  
  # Long pivot
  mapping_expanded <- mapping %>%
    mutate(ortholog = strsplit(ortholog, ";")) %>%
    tidyr::unnest(ortholog) %>%
    mutate(ortholog = trimws(ortholog))
  
  # Logical array
  indices <- which(mapping_expanded$daphnia %in% genes)
  
  
  return(mapping_expanded$ortholog[indices])
}


# Peak to Metabolte
mz_to_metabolite <- function(mzs, mapping) {
  
  # Unique genes
  mzs <- unique(mzs)
  
  # Long pivot
  mapping_expanded <- mapping %>%
    mutate(metabolite = strsplit(metabolite, ";")) %>%
    tidyr::unnest(metabolite) %>%
    mutate(metabolite = trimws(metabolite))
  
  # Logical array
  indices <- which(mapping_expanded$mz %in% mzs)
  
  
  return(mapping_expanded$metabolite[indices])
}

#
#