fit.m1 = function(udata,temp.info,cand.tref=55:80,dur=5){
  ## dur is duration. if we exclude weekends dur=5 means one week
  ## udata is n by 24 usage data. assume it is in order. for ex, 1st row is Jan-01-2011, 2nd row is Jan-02-2011
  ## temp.info is n by 24 temperature data. assume it is aligned with udata
  ## cand.tref is a list of candidates of Tr
  ## m1 is U(t, d) = A+(t) * (T(t, d) - Tr(w(d)))+  +  A-(t) * (Tr(w(d)) - T(t,d))+  + C(t) + e
  
  n = nrow(udata)
  ## at first, set Tr(w(d)) among 55~80
  ## start from U(t, d) = A+(t) * (T(t, d) - Tr)+  +  A-(t) * (Tr - T(t,d))+  + C(t) + e
  tmp.rss = c()
  for (i in cand.tref){
    res = try.tref(udata,temp.info,i) 
    tmp.rss = c(tmp.rss,res$rss)
  }
  tref = cand.tref[which.min(tmp.rss)] ## now set coefs starting point
  Trseq = rep(tref,ceiling(n/dur))
  estep = try.tref(udata,temp.info,tref) 
  
  mstep = fit.tref.seq(estep$coefs,udata,temp.info,cand.tref,dur)
  estep = fit.coefs(mstep$Trseq,udata,temp.info,dur)
  prev.mstep = mstep
  prev.estep = estep
    
  ## start while loop for EM algorithm
  while(1){
    ## with fixed coefs find Tr(w(d)) sequence
    mstep = fit.tref.seq(estep$coefs,udata,temp.info,cand.tref,dur)
    estep = fit.coefs(mstep$Trseq,udata,temp.info,dur)
    if (prev.estep$rss<=estep$rss) break
    prev.mstep = mstep
    prev.estep = estep    
  }  
  list(rss=prev.estep$rss,est=prev.estep$est,coefs=prev.estep$coefs, trseq=prev.mstep$Trseq)
}

## for each period, find Tref
fit.tref.seq=function(coefs,udata,temp.info,cand.tref,d){
  cand.ref.t = cand.tref
  coefs[is.na(coefs)]=0
  cnt=0
  datalen = nrow(udata)
  trseq = matrix(0,1,ceiling(datalen/d))
  rss = matrix(0,1,ceiling(datalen/d))
  for (j in 1:ceiling(datalen/d)){
    ls = min(d,datalen-d*j+d); err = c()
    for (i in cand.ref.t){
      tmp.err = 0
      for (k in 1:ls){
        tmp.err = tmp.err + sum((udata[d*cnt+k,]-coefs[seq(2,72,3)]*pmax(temp.info[d*cnt+k,]-i,0)-coefs[seq(3,72,3)]*pmax(i-temp.info[d*cnt+k,],0)-coefs[seq(1,72,3)])^2,na.rm=T)
      } 
      err = c(err,tmp.err)
    }
    cnt = cnt+1
    trseq[cnt] = cand.ref.t[which.min(err)]
    rss[cnt] = min(err)
    cand.ref.t = (trseq[cnt]-2):(trseq[cnt]+2) ## plus minus 2 allowed
    cand.ref.t = cand.ref.t[cand.ref.t>=min(cand.tref) & cand.ref.t<=max(cand.tref)]
  }
  list(Trseq = trseq, rss = sum(rss))
}

## given Tref sequence, find coefs
fit.coefs=function(ref.temp,udata,temp.info,d){
  fit.data= c();cnt=0
  datalen = nrow(udata)
  for (j in 1:ceiling(datalen/d)){
    ls = min(d,datalen-d*j+d)
    for (k in 1:ls){
      cnt = cnt+1
      dummy1 = matrix(0,24,24);dummy2 = matrix(0,24,24)
      diag(dummy1)= pmax(temp.info[cnt,]-ref.temp[j],0)
      diag(dummy2)= pmax(ref.temp[j]-temp.info[cnt,],0)
      fit.data = rbind(fit.data, cbind(dummy1,dummy2,t(udata[cnt,])))
    }
  }
  coefs=matrix(0,24,3);rss=0;est=matrix(0,datalen,24)
  for (k in 1:24){
    fit2data = fit.data[seq(k,24*datalen,24),c(k,k+24,49)]
    fit2data = data.frame(fit2data); names(fit2data)[3]='y'
    fit = lm(y~.,data=fit2data)
    est[,k] = fit$fit
    coefs[k,] = fit$coef
    rss = rss + sum(fit$res^2,na.rm=T)
  }
  coefs[is.na(coefs)]=0
  list(coefs=coefs, rss=rss, est=est)
}

try.tref=function(udata,temp.info,tref){
  ## fit U(t, d) = A+(t) * (T(t, d) - Tr)+  +  A-(t) * (Tr - T(t,d))+  + C(t) + e with Tr = tref
  fit.data = c();n = nrow(udata) 
  for (j in 1:n){
    dummy1 = matrix(0,24,24); dummy2 = matrix(0,24,24)
    diag(dummy1)=pmax(temp.info[j,]-tref,0); diag(dummy2)=pmax(tref-temp.info[j,],0)
    fit.data = rbind(fit.data, cbind(dummy1,dummy2,t(udata[j,])))
  }
  rss = 0;est=c();coefs=c()
  for (k in 1:24){
    fit2data = fit.data[seq(k,24*n,24),c(k,k+24,49)]
    fit2data = data.frame(fit2data); names(fit2data)[3]='y'
    fit = lm(y~.,data=fit2data)
    est = rbind(est,fit$fit)
    coefs=c(coefs,fit$coef)
    rss = rss + sum(fit$res^2,na.rm=T)
  }
  coefs[is.na(coefs)]=0
  list(rss=rss,coefs=coefs,est=est)
}  

fit.m2 = function(udata,temp.info,cand.tref=55:80){
  ## udata is n by 24 usage data. assume it is in order. for ex, 1st row is Jan-01-2011, 2nd row is Jan-02-2011
  ## temp.info is n by 24 temperature data. assume it is aligned with udata
  ## cand.tref is a list of candidates of Tr
  ## m2 is U(t, d) = A+(t) * (T(t, d) - Tr(t))+  +  A-(t) * (Tr(t) - T(t,d))+  + C(t) + e
  ## m2 is time wise model. Just fit 24 times
  n = nrow(udata)
  
  coefs = matrix(0,24,3)
  tmp.rss = matrix(0,1,24)
  rss = matrix(10000,1,24)
  est = matrix(0,n,24)
  tref = matrix(0,1,24)
  for (i in cand.tref){
    fit.data = c()
    for (j in 1:n){
      dummy1 = matrix(0,24,24); dummy2 = matrix(0,24,24)
      diag(dummy1)=pmax(temp.info[j,]-i,0); diag(dummy2)=pmax(i-temp.info[j,],0)
      fit.data = rbind(fit.data, cbind(dummy1,dummy2,t(udata[j,])))
    }
    for (k in 1:24){
      fit2data = fit.data[seq(k,24*n,24),c(k,k+24,49)]
      fit2data = data.frame(fit2data); names(fit2data)[3]='y'
      fit = lm(y~.,data=fit2data)
      tmp.rss = sum(fit$res^2,na.rm=T)
      if (tmp.rss < rss[k]) {
        rss[k] = tmp.rss
        coefs[k,] = fit$coef
        est[,k] = fit$fit
        tref[k] = i
      }
    }
  }
  coefs[is.na(coefs)]=0
  list(rss=sum(rss),est=est,coefs=coefs,trseq=tref)
}

fit.m3 = function(udata,temp.info,cand.tref=55:80,d=5,coefs){
  ## d is duration. if we exclude weekends d=5 means one week
  ## udata is n by 24 usage data. assume it is in order. for ex, 1st row is Jan-01-2011, 2nd row is Jan-02-2011
  ## temp.info is n by 24 temperature data. assume it is aligned with udata
  ## cand.tref is a list of candidates of Tr
  ## coefs is from m2 result. 24 by 3 matrix
  ## m3 is U(t, d) = A+(t) * (T(t, d) - Tr(w(d),t))+  +  A-(t) * (Tr(w(d),t) - T(t,d))+  + C(t) + e
  ## m3 starts from m2. fix coefs first for each t, and find Tr seq for each t
  ## so, run m2, m3  sequentially

  n = nrow(udata)
  prev.rss = 10000
  trseq = matrix(0,24,ceiling(n/d))
  while(1){
    rss=0
    for (k in 1:24){
      cand.ref.t = cand.tref
      cnt=0

      for (j in 1:ceiling(n/d)){
        ls = min(d,n-d*j+d); err = c()
        for (i in cand.ref.t){
          err = c(err,sum((udata[d*cnt+1:ls,k]-coefs[k,2]*pmax(temp.info[d*cnt+1:ls,k]-i,0)-coefs[k,3]*pmax(i-temp.info[d*cnt+1:ls,k],0)-coefs[k,1])^2,na.rm=T))
        }
        cnt = cnt+1
        trseq[k,j] = cand.ref.t[which.min(err)]
        rss = rss + min(err)
        cand.ref.t = (trseq[k,j]-2):(trseq[k,j]+2) ## plus minus 2 allowed
        cand.ref.t = cand.ref.t[cand.ref.t>=min(cand.tref) & cand.ref.t<=max(cand.tref)]
      }
    }
    ## trseq is set now
  
    ## need to fit coefs with fixed trseq
    fit.data= c();cnt=0
    for (j in 1:ceiling(n/d)){
      ls = min(d,n-d*j+d)
      for (k in 1:ls){
        cnt = cnt+1
        dummy1 = matrix(0,24,24);dummy2 = matrix(0,24,24)
        diag(dummy1)= pmax(temp.info[cnt,]-trseq[,j],0)
        diag(dummy2)= pmax(trseq[,j]-temp.info[cnt,],0)
        fit.data = rbind(fit.data, cbind(dummy1,dummy2,t(udata[cnt,])))
      }
    }
    coefs=matrix(0,24,3);rss=0;est=c()
    for (k in 1:24){
      fit2data = fit.data[seq(k,24*n,24),c(k,k+24,49)]
      fit2data = data.frame(fit2data); names(fit2data)[3]='y'
      fit = lm(y~.,data=fit2data)
      est = rbind(est,fit$fit)
      coefs[k,]=fit$coef
      rss = rss + sum(fit$res^2,na.rm=T)
    }
    if (prev.rss <= rss) break
    prev.rss = rss
    prev.est = est
    prev.coefs = coefs
    prev.trseq = trseq
  }
  prev.coefs[is.na(prev.coefs)]=0
  list(rss=prev.rss,est=t(prev.est),coefs=prev.coefs,trseq=prev.trseq)
}
