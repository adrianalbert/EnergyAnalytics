ksc.dist = function(x,y){
  if (sum(abs(y))!=0 & sum(abs(x)!=0)) {
    alpha = sum(x*y)/sum(y^2)
    return(sum((x-alpha*y)^2)/sum(x^2))
  } else return(1)
}

dhat.shift = function(x,y){
  min.d = ksc.dist(x,y)
  
  len = length(y)
  range = 0
  opt.shift = 0; opt.y = y
  for (shift in range){
    y.shift = y
    if (shift<0) y.shift = c(y[(1-shift):len], rep(0,-shift))
    else if (shift>0) y.shift = c(rep(0,shift), y[1:(len-shift)])
    
    cur.d = ksc.dist(x,y.shift)
    if (cur.d < min.d){
      opt.shift = shift; opt.y = y.shift; min.d = cur.d
    }    
  }
  list(min.d = min.d, opt.shift = opt.shift, opt.y = opt.y)
}

ksc.center.update = function(mem,A,k,cur.center=NULL){
  ## mem: n by 1, Membership for each time series.
  ## A: n by p, each row is a time series
  ## cur.center: k by p, current cluster centeroids
  n = dim(A)[1]; p = dim(A)[2]
  opt.A = c(); new.center = c()
  
  if (is.null(cur.center)){
    opt.A = A
  } else {
    for (i in 1:n){
      opt.A = rbind(opt.A, dhat.shift(A[i,],cur.center[mem[i],])$opt.y)
    }
  }
  for (i in 1:k){
    if (sum(mem==i)==0){
      new.center = rbind(new.center, rep(0,p))
    } else {
      tmp.A = matrix(opt.A[mem==i,],ncol=p)
      b = tmp.A / matrix(rep(sqrt(apply(tmp.A^2,1,sum)),p),ncol=p)
      M = t(b)%*%b - dim(tmp.A)[1]*diag(p)
      ks = eigen(M)$vectors[,1]     
      
      if (sum(ks)<0) ks = -ks
      new.center = rbind(new.center, ks)
    }
  }
  new.center
}

ksc.assign.update = function(A,k,cur.center){
  n = dim(A)[1]; mem=c()
  for (i in 1:n){
    min.d = sum(A^2); min.mem = 0
    for (j in 1:k){
      tmp.d = dhat.shift(A[i,],cur.center[j,])$min.d
      if (tmp.d < min.d) {
        min.d = tmp.d
        min.mem = j
      }
    }
    mem[i] = min.mem
  }
  mem
}

ksc = function(A,k,init.cen=NULL,init.mem=NULL,max.iter=100){
  n = dim(A)[1]; p = dim(A)[2]
  ## centeroids init
  if (!is.null(init.cen)) {
    cur.cen = init.cen
  } else if (!is.null(init.mem)) {
    cur.cen = ksc.center.update(init.mem,A,k)
  } else {cur.cen = A[sample(n,k),]}
  
  ## assignment init
  cur.mem = ksc.assign.update(A,k,cur.cen)
  
  for (i in 1:max.iter){
    prev.mem = cur.mem
    cur.cen = ksc.center.update(cur.mem, A,k)
    cur.mem = ksc.assign.update(A,k,cur.cen)
    if (sum(abs(prev.mem-cur.mem))==0) break
  }
  list(mem = cur.mem, center = cur.cen)
}

house.dwt = function(x){
  ## assume length of x is 2^l
  len = length(x)
  prev.d = x
  coefs = c()
  while(len>1){
    tmp.d = (prev.d[seq(1,len,2)]-prev.d[seq(2,len,2)])/2
    coefs = c(tmp.d, coefs)
    prev.d = (prev.d[seq(1,len,2)]+prev.d[seq(2,len,2)])/2
    len = len/2
  }
  c(prev.d, coefs)
}

house.idwt = function(coefs,level){
  ## assume length of coefs is 2^l
  ## level: level of reconstruction: level l will return 2^l data points.
  ## and level should be equal or less than l
  len = length(coefs)
  x = coefs[1]
  for (i in 1:level) x = matrix(rbind(x+coefs[(2^(i-1)+1):2^i], x-coefs[(2^(i-1)+1):2^i]),nrow=1)
  x
}

inc.ksc = function(A,k,start.level=4,max.ter=100){
  ## incremental ksc will be used in case p is very big
  ## assume that p is 2^l
  n = dim(A)[1]; p = dim(A)[2]; y = c()
  z = t(apply(A,1,house.dwt))
  cur.mem = NULL
  for (i in start.level:log2(p)){
    y = t(apply(z,1,house.idwt,i))
    res = ksc(A,k,init.mem = cur.mem)
    cur.mem = res$mem
    cur.cen = res$center
  }
  list(mem = cur.mem, center = cur.cen)
}