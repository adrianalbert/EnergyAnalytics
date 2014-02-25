# clustering_wrappers.r
#
# Class to encapsulate clustering analysis on occupancy states. 
# 
# Adrian Albert
# Last modified: April 2013.
# -----------------------------------------------------------------------

# ________________________________
# Hard K-Means clustering wrapper

kmeans_wrapper = function(data, Kmin = 3, Kmax = 3) {
  # perform clustering
  cat('k = ')
  if (Kmax > Kmin) {
    obj  = rep(0, Kmax-Kmin+1)
    expl = rep(0, Kmax-Kmin+1)
    Kvec = Kmin:Kmax
    for (k in 1:length(Kvec)) {
      if (verbose) cat(paste(Kvec[k], '; '))
      fit        = kmeans(data, Kvec[k], nstart = 10, iter.max = 100)
      expl[k]    = fit$betweenss / fit$totss
      obj[k]     = fit$tot.withinss
    }            
    if (verbose) cat('\n')
    names(obj)   = as.character(Kmin:Kmax)
    names(expl)  = as.character(Kmin:Kmax)
    return(list(Objective = obj, Explained = expl))
  } else {
    fit          = kmeans(data, Kmin, nstart = 10, iter.max = 100)       
    return(list(centers = fit$centers, assign = fit$cluster))
  }            
}

# _________________________________
# Soft K-Means clustering wrapper

library('e1071')
cmeans_wrapper = function(data, Kmin = 3, Kmax = 3) {
  # perform clustering
  cat('k = ')
  if (Kmax > Kmin) {
    obj  = rep(0, Kmax-Kmin+1)
    Kvec = Kmin:Kmax
    for (k in 1:length(Kvec)) {
      if (verbose) cat(paste(Kvec[k], '; '))
      fit        = cmeans(data, Kvec[k], iter.max = 100)
      obj[k]     = fit$withinerr
    }            
    if (verbose) cat('\n')
    names(obj)   = as.character(Kmin:Kmax)
    return(obj)
  } else {
    fit          = cmeans(data, Kmin, iter.max = 100)
    return(list(centers = fit$centers, assign = fit$cluster, membership = fit$membership))
  }            
}

