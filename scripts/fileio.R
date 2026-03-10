load.dataset <- function(meta.file, meta.sep = ',', data.file, data.sep = ','){
  meta <- read.table(meta.file, header = T, sep = meta.sep, row.names = 1, check.names = F)
  data <- read.table(data.file, header = T, sep = data.sep, row.names = 1, check.names = F)

  data <- as.matrix(data)

  mids <- match(row.names(meta), colnames(data))
  data <- t(data[,mids[!is.na(mids)]])
  meta <- meta[!is.na(mids),]
  
  stopifnot(row.names(meta) == row.names(data))
  stopifnot(!any(is.na(row.names(meta))))
  stopifnot(!any(is.na(row.names(data))))
  
  return(list(meta.data = meta, data.matrix = data))
}
