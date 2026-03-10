rm(list = ls())


library(WGCNA)
source('fileio.R')


# allow multi-threading within WGCNA
# skip this line if you run RStudio or other third-party R environments
#enableWGCNAThreads()


# load data
ds = load.dataset(
  meta.file = '../data/sample_sheet.csv', meta.sep = ',',
  data.file = '../data/rna_norm_counts.csv', data.sep = ','
)
data = ds$data.matrix
dim(data)


# search soft-thresholding powers
powers = 2:20
sft = pickSoftThreshold(data, powerVector = powers, verbose = 5)

# plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = 'Soft Threshold (power)', ylab = 'Scale Free Topology Model Fit,signed R^2', type = 'n',
     main = paste('Scale independence'));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = 0.9, col = 'red');
abline(h = 0.90, col = 'red')
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = 'Soft Threshold (power)', ylab = 'Mean Connectivity', type = 'n',
     main = paste('Mean connectivity'))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.9, col = 'red')


# apply soft-thresholding
softPower = 4
adjacency = adjacency(data, power = softPower)


# create topological overlap matrix
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM


# hierarchical clustering for community detection
geneTree = hclust(as.dist(dissTOM), method = 'average')
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 20)
moduleColors = labels2colors(dynamicMods)
summary(as.factor(moduleColors))

# plot the results
sizeGrWindow(12,9)
plotDendroAndColors(geneTree, moduleColors, 'Dynamic Tree Cut',
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


# calculate eigengenes
MEList = moduleEigengenes(data, colors = moduleColors)
MEs = MEList$eigengenes

# cluster module eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");

sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")


# correlate eigengenes with labels
conds <- lapply(ds$meta.data, function(x) as.numeric(as.factor(x)))
conds = as.data.frame(conds, row.names = row.names(ds$meta.data))

moduleTraitCor = cor(MEs, conds, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(data))

# plot the results
sizeGrWindow(10,6)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(conds),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


# export to cytoscape 
modules = c('cyan', 'tan') # select to modules of interest
probes = colnames(data)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
modTOM = 1 - dissTOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
efname = paste('../data/', paste(modules, collapse = '-'), '_edges.txt', sep = '')
nfname = paste('../data/', paste(modules, collapse = '-'), '_nodes.txt', sep = '')

exportNetworkToCytoscape(modTOM,
                         edgeFile = efname,
                         nodeFile = nfname,
                         weighted = TRUE, threshold = 0.25, 
                         nodeNames = modProbes,
                         nodeAttr = moduleColors[inModule])
