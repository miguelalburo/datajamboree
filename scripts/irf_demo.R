rm(list = ls())

library(iRF)
library(AUC)
source('fileio.R')


# load data
ds <- load.dataset(
  meta.file = '../data/sample_sheet.csv', meta.sep = ',',
  data.file = '../data/rna_norm_counts.csv', data.sep = ','
)
X <- as.data.frame(ds$data.matrix)
Y <- ds$meta.data$REF


# prepare inputs
X <- X[Y != 'Control',]
Y <- as.factor(Y[Y != 'Control'])
Y <- as.factor(as.numeric(relevel(Y, '1x'))-1)

n <- dim(X)[1]
p <- dim(X)[2]


# split training-testing set
train.id <- sample(1:n, size = 8*round(n/10))
test.id <- setdiff(1:n, train.id)


# fit iRF without interaction importance estimation
sel.prob <- rep(1/p, p)

fit <- iRF(x = X[train.id,], 
           y = Y[train.id], 
           xtest = X[test.id,], 
           ytest = Y[test.id],
           n.iter = 5, 
           iter.return = 1:5,
           n.core = 4
)


# plot ROC/AUC
rf <- fit$rf.list

plot(0:1, 0:1, type = 'l', lty = 2, xlab = 'FPR', ylab = 'TPR', main = 'ROC Curve')
for (iter in 1:5){
  cat(paste('iter = ', iter, ':: '))
  roc.info <- roc(rf[[iter]]$test$votes[,2], Y[test.id])
  lines(roc.info$fpr, roc.info$tpr, type = 'l', col = iter, lwd = 2)
  cat(paste('AUROC: ', round(100*auc(roc.info), 2), '%\n', sep = ''))
}
legend('bottomright', legend=paste('iter:', 1:iter), col=1:iter, lwd=2, bty='n')


# plot feature weights
par(mfrow = c(1,5))
for (iter in 1:5){
  varImpPlot(rf[[iter]], n.var = 10, main = paste('Variable Importance (iter:', iter, ')')) 
}


# fit iRF with interaction importance estimation
fit <- iRF(x = X[train.id,], 
           y = Y[train.id], 
           xtest = X[test.id,], 
           ytest = Y[test.id],
           n.iter = 5, 
           iter.return = 1:5,
           select.iter = T,
           # int.return = 2:5,
           n.bootstrap = 50, 
           n.core = 4
)


# plot iteration stability scores
rf.ints <- fit$interaction 
# rf.ints <- fit$interaction[[2]]

toplot <- rf.ints$stability
names(toplot) <- rf.ints$int
toplot <- sort(toplot, decreasing = T)

dev.off()
dotchart(rev(toplot[1:min(20, length(toplot))]), 
         xlab = 'Stability Score', 
         xlim = c(0, 1)
)
