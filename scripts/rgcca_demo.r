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


# load DIMS data
polar.ds <- load.dataset(
    meta.file = '../data/sample_sheet.csv', meta.sep = ',',
    data.file = '../data/polar_pos_pqn_imputed_glog.csv', data.sep = ','
)


# take samples available on both data
shared.ids <- intersect(row.names(rna.ds$data.matrix), row.names(polar.ds$data.matrix))
A <- list(rna.ds$data.matrix[shared.ids,], polar.ds$data.matrix[shared.ids,])


# apply RGCCA/SGCCA
A <- lapply(A, function(x) scale(x))
C <- matrix(c(0, 1, 1, 0), 2, 2)

cca.res <- rgcca(A, C, tau = c(1,1), ncomp = c(2,2), scale = FALSE, verbose = FALSE)
cca.res$AVE

# if you are using an older version RGCCA
cca.res <- sgcca(A, C, c1 = c(0.2,0.2), ncomp = c(2,2), scale = FALSE, verbose = FALSE)
# if you are using a newer version RGCCA
cca.res <- rgcca(A, C, sparsity = c(0.2,0.2), ncomp = c(2,2), scale = FALSE, verbose = FALSE)
cca.res$AVE


# number of features selected
colSums(cca.res$a[[1]] != 0)
colSums(cca.res$a[[2]] != 0)


# CCA plots
factor.plot <- function(y1, y2, conds){
  df <- data.frame(datablock1 = y1,  datablock2 = y2, colFactor = conds)
  p <- ggplot(df, aes(datablock1, datablock2)) + 
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
    ggtitle('Factor plot') + 
    geom_text(aes(colour = colFactor, label = rownames(df)), vjust = 0, nudge_y = 0.03, size = 3) + 
    theme(legend.position = 'bottom', legend.box = 'horizontal', legend.title = element_blank())
  return(p)
}

group <- rna.ds$meta.data[shared.ids,]$REF
p <- factor.plot(cca.res$Y[[1]][,1], cca.res$Y[[2]][,1], group)
ggplotly(p)

group <- rna.ds$meta.data[shared.ids,]$Site
p <- factor.plot(cca.res$Y[[1]][,2], cca.res$Y[[2]][,2], group)
ggplotly(p)
