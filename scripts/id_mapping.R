rm(list = ls())


# read target IDs
# df <- read.csv("../data/rna_norm_counts.csv", header = T, row.names = 1, check.names = F)
df <- read.csv('../data/polar_pos_pqn_imputed_glog.csv', header = T, row.names = 1, check.names = F)

ids <- row.names(df)
length(ids)


# load mapping files
# id.maps <- read.csv("../annotations/rna_dma_to_dpx_mappings.tsv", header = F, row.names = 1, sep = '\t', stringsAsFactors = F)
id.maps <- read.csv("../annotations/polar_pos_pkl_to_kegg_annotations.tsv", header = F, row.names = 1, sep = '\t', stringsAsFactors = F)

dim(id.maps)
head(id.maps)


# mapping
id.maps.matched <- subset(id.maps, row.names(id.maps) %in% ids)
id.maps.matched <- id.maps.matched[id.maps.matched$V3 != '',]
dim(id.maps.matched)

mapped.ids <- unique(unlist(strsplit(id.maps.matched$V3, split = ';')))
length(mapped.ids)


# save to file
writeLines(mapped.ids, '../data/polar_pos_pqn_imputed_glog_kegg_ids.csv')

