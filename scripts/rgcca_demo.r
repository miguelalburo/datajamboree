rm(list = ls())


library(RGCCA)
library(ggplot2)
library(plotly)
source('fileio.R')


# load RNA-seq data
rna.ds <- load.dataset(
    meta.file = '../data/sample_sheet.csv', meta.sep = ',',
    data.file = '../data/rna_norm_counts.csv', data.sep = ','
)


# load DIMS data | positive ion mode 149 samples
polar.ds <- load.dataset(
    meta.file = '../data/sample_sheet.csv', meta.sep = ',',
    data.file = '../data/polar_pos_pqn_imputed_glog.csv', data.sep = ','
)


# take samples available on both data
shared.ids <- intersect(row.names(rna.ds$data.matrix), row.names(polar.ds$data.matrix))
# list of peak intensity / gene counts matrices
A <- list(rna.ds$data.matrix[shared.ids,], polar.ds$data.matrix[shared.ids,])


# apply RGCCA/SGCCA
A <- lapply(A, function(x) scale(x))
C <- matrix(c(0, 1, 1, 0), 2, 2)

cca.res <- rgcca(A, C, sparsity = c(0.2,0.2), ncomp = c(2,2), scale = FALSE, verbose = FALSE)


cca.res$AVE 
# Average Variance Explained
# Block 1 = RNA seq data, Block 2 = DIMS data
# AVE inner = correlation between the components of the two blocks
# AVE outer = correlation between the original variables and the components of each block
# Output: Component 1 | Component 2



# number of features selected 
# $a gives a list. Each element is a matrix of block weights
colSums(cca.res$a[[1]] != 0) # all genes 18903
colSums(cca.res$a[[2]] != 0) # all metabolites 1285


# CCA plots

factor.plot <- function(y1, y2, conds){
  df <- data.frame(datablock1 = y1,  datablock2 = y2, color = conds)
  p <- ggplot(df, aes(datablock1, datablock2, color = color)) + 
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
    ggtitle('Factor plot') + 
    geom_point(aes(size = 2)) + 
    theme(legend.position = 'bottom', legend.box = 'horizontal', legend.title = element_blank()) +
    theme_bw()
  return(p)
}

# cca.res$Y[[1]][,1] : block components, block 1 component 1

group <- rna.ds$meta.data[shared.ids,]$REF
p <- factor.plot(cca.res$Y[[1]][,1], cca.res$Y[[2]][,1], group)
print(p)

group <- rna.ds$meta.data[shared.ids,]$REF
p <- factor.plot(cca.res$Y[[1]][,2], cca.res$Y[[2]][,2], group)
print(p)

# cca.res$Y[[1]][,2] : block components, component 2
group <- rna.ds$meta.data[shared.ids,]$Site
p <- factor.plot(cca.res$Y[[1]][,2], cca.res$Y[[2]][,2], group)
print(p)

group <- rna.ds$meta.data[shared.ids,]$Site
p <- factor.plot(cca.res$Y[[1]][,1], cca.res$Y[[2]][,1], group)
print(p)


