#
# Regularised Generalised Canonical Correlation Analysis (RGCCA)
#
# Design
#
# DIMS Peak Intensity ~ RNAseq Gene Counts ~ Chemical Concentrations
#


# Setup -------------------------------------------------------------------

# Essentials
source('functions.R')


# Specifics
library(RGCCA)
library(plotly)
source('fileio.R')


# Load data

## Gene Counts | Normalised
rna.ds <- load.dataset(
    meta.file = '../data/sample_sheet.csv', meta.sep = ',',
    data.file = '../data/rna_norm_counts.csv', data.sep = ','
)


## DIMS Peak Intensity | Positive Ion Mode
polar.ds <- load.dataset(
    meta.file = '../data/sample_sheet.csv', meta.sep = ',',
    data.file = '../data/polar_pos_pqn_imputed_glog.csv', data.sep = ','
)


# Combining into Feature Matrix List
shared.ids <- intersect(row.names(rna.ds$data.matrix), row.names(polar.ds$data.matrix))
A <- list(rna.ds$data.matrix[shared.ids,], polar.ds$data.matrix[shared.ids,])


# Analysis ----------------------------------------------------------


# Scaling
A <- lapply(A, function(x) scale(x))


# Design Matrix
C <- matrix(c(0, 1, 1, 0), 2, 2)


# Perform RGCCA
sparsity = c(0.2,0.2)
tau = "optimal"
ncomp = 2 # no. components to calculate

cca.res <- rgcca(
  blocks = A, 
  connection = C,
  sparsity = sparsity,
  tau = tau,
  ncomp = ncomp,
  method = "rgcca",
  scale = FALSE, verbose = FALSE
)


# Results -----------------------------------------------------------------


# Average Variance Explained
#
# Block 1 = RNA seq data, Block 2 = DIMS data
# AVE inner = between blocks
# AVE outer = within blocks
# Output: Component 1 | Component 2

cca.res$AVE


# Number of features selected 
#
# $a gives a list. Each element is a matrix of block weights
# Output: Component 1 | Component 2

colSums(cca.res$a[[1]] != 0) # 1324 1234
colSums(cca.res$a[[2]] != 0) # 91 81


# List of genes which are selected in the either components
# UNION approach
genes_of_interest <- cca.res$a[[1]] %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>% 
  filter(V1 != 0 | V2 != 0) %>%
  pull(gene)

# List of metabolites which are selected in the either components
# UNION approach
metabolites_of_interest <- cca.res$a[[2]] %>%
  as.data.frame() %>%
  rownames_to_column('metabolite') %>% 
  filter(V1 != 0 | V2 != 0) %>%
  pull(metabolite)


# Plots -------------------------------------------------------------------


# Component Factor Plot 
conds <- c("Site", "REF", "Description")
p <- list()

for (cond in conds){
  for (component in seq(ncomp)){
    
    group <- rna.ds$meta.data[shared.ids,][[cond]]
    y1 <- cca.res$Y[[1]][,component]
    y2 <- cca.res$Y[[2]][,component]
    
    name <- paste(cond, "Component", component, sep = "_")
    p[[name]] <- factor.plot(y1, y2, group) + 
      labs(subtitle = paste("Component", component))
  }
}



# Saving & Writing --------------------------------------------------------


# Write genes of interest to .txt
write.table(genes_of_interest, '../results/rgcca/gois_rgcca.txt',
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Write metabolites of interest to .txt
write.table(metabolites_of_interest, '../results/rgcca/mois_rgcca.txt',
            row.names = FALSE, col.names = FALSE, quote = FALSE)


# Writing plots to .pdf
for (fplot in names(p)){
  pdf(paste('../results/rgcca/factor_plots_', fplot, '.pdf', sep = ''), width = 10, height = 6)
  print(p[[fplot]])
  dev.off()
}