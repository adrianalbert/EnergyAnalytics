library('gbm')
gbmImpute = function(x, max.iters = 2, cv.fold = 5, verbose=T, n.trees = 200) {
  
  missing.matrix = is.na(x)
  numMissing = sum(missing.matrix)
  if(verbose) {
    print(paste("imputing on", numMissing, "missing values with matrix size",
                nrow(x)*ncol(x), sep=" "))
  }
  if(numMissing == 0) {
    return (x)
  }
  
  missing.cols.indices = which(apply(missing.matrix, 2, function(i) {
    any(i)
  }))
  
  for (i in 1:max.iters) {
    if (verbose) print (paste("Begin iteration: ", i))
    x[,missing.cols.indices] = sapply(missing.cols.indices, function(j) {
      good.data = which(!missing.matrix[,j])
      gbm1 <- gbm(x[good.data,j] ~ .,
                  data = as.data.frame(x[good.data,-j]),
                  var.monotone = rep(0, ncol(x)-1), # -1: monotone decrease,
                  # +1: monotone increase,
                  #  0: no monotone restrictions
                  distribution="gaussian",     # bernoulli, adaboost, gaussian,
                  # poisson, coxph, and quantile available
                  n.trees=n.trees,                # number of trees
                  shrinkage=0.005,             # shrinkage or learning rate,
                  # 0.001 to 0.1 usually work
                  interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc.
                  bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                  train.fraction = 0.5,        # fraction of data for training,
                  # first train.fraction*N used for training
                  n.minobsinnode = 10,         # minimum total weight needed in each node
                  cv.folds = cv.fold,                # do 5-fold cross-validation
                  keep.data=TRUE,              # keep a copy of the dataset with the object
                  verbose=verbose)                # print out progress
      best.iter <- gbm.perf(gbm1,method="OOB", plot.it = F)
      data.predict = predict(gbm1, newdata = as.data.frame(x[-good.data,-j]), n.trees = best.iter)
      x[-good.data,j] = data.predict
      x[,j]
    })
  }
  
  return ( list (
    x=x,
    missing.matrix=missing.matrix
  ))
}
