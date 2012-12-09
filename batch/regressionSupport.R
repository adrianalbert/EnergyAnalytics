# returns a matrix of length(sort(unique(membership))) columns, that is all 
# zeros except for the regressor values that match the membership values 
# corresponding to the column number this supports regression with separate 
# coefficints for each sub group defined by the membership
# for example, splitRegressor(Tout,dates$hour) will return a matrix with 24 
# columns, where the only non-zero entry per row contains the Tout value in
# the column corresponding to the hour of day it was recorded 
regressor.split = function(regressor,membership) {
  mat <- c()
  nm  <- c()
  for (i in sort(unique(membership))) {
    mat <- cbind(mat,ifelse(membership==i,1,0)) # add a colunm of 1's and 0's
    nm <- c(nm,i)
  }
  colnames(mat) <- nm
  return(mat * regressor)
}

# break a vector out for continuous piecewise regression (i.e. the fitted
# segments will join eachother at each end) into a matrix whose row
# totals are the original values, but whose columns divide the value across
# a set of bins, so 82 across bins with boundaries c(50,60,65,70,80,90) 
# becomes the row 50,10,5,5,10,2,0, which sums to 82...
# This is very useful for finding rough change points in thermal response
# as is expected for buildings with clear setpoints

# TODO: This can create a column of zeros, which should break the regression
# so we might need to prune the columns when we're done and keep track of
# which bins are in play when comparing across regressions
regressor.piecewise = function(regressor,bins) {
  binLower = 0
  mat <- c()
  nm = c()
  for(binUpper in c(bins,Inf)) {
    col = regressor - binLower
    col[col < 0] = 0
    col[(col> binUpper-binLower)] = binUpper-binLower
    mat = cbind(mat,col)
    nm = c(nm,paste(binLower,'-',binUpper,sep=''))
    binLower = binUpper
  }
  colnames(mat) <- nm
  return(mat)
}