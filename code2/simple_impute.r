## imputation function
sim.impute = function(A,uidx=4:99){
  require('imputation')
## assume the input comes for the same user and ordered by date
## at least two rows should be valid
## uidx = 4:99 is column indexes for usage data

  n = nrow(A);need = 0 ## need to run knn impute
  im.idx = which(apply(A[,uidx],1,function(i){
    if (sum(is.na(i))>0) return(1)
    else if (sum(i)==0) return(1)
    else return(0)
  })==1)
  
  ## return fail if all values are zero
  if (sum(A[,uidx],na.rm=T)==0) {print("all values are zero=>can't impute");return(A)}
  
  for (i in im.idx){
    if (any(is.na(A[i,uidx]))){
      ## if it's long: all null
      if (sum(is.na(A[i,uidx]))==length(uidx)){
        ## put the mean or mean of 2 near days
        ## A[i,uidx] = apply(A[-im.idx,uidx],2,mean) ## mean of all data, not good
        near = order(abs(1:n-i))
        A[i,uidx] = apply(A[near[!(near%in%im.idx)][1:2],uidx],2,mean)
      } else {need = 1}
      ## if null value length is short
      ## leave them and use imputation package.
    } else if (sum(A[i,uidx])==0){ ## all zero rows
      near = order(abs(1:n-i))
      A[i,uidx] = apply(A[near[!(near%in%im.idx)][1:2],uidx],2,mean)
    } 
  }
  if (need) A[,uidx] = kNNImpute(A[,uidx], k=3, verbose=F)$x
  A
}
