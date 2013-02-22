#!/usr/bin/Rscript

#.libPaths('~/R/library') # use my local R library even from the command line

require(ggplot2)
require(gtools)
require(gridExtra)
require(reshape2)
require(plyr)

# todo: is there a better way to detect the current directory?
conf.basePath = file.path('~/EnergyAnalytics/batch')
if(Sys.info()['sysname'] == 'Windows') {
  conf.basePath = file.path('f:/dev/pge_collab/EnergyAnalytics/batch')
}
setwd(conf.basePath)

resultsDir = 'results_climate_test' # real
#resultsDir = 'test_results' # sub -sample for testing

# run 'source' on all includes to load them 
source(file.path(getwd(),'localConf.R'))         # Sam's local computer specific configuration
source(file.path(getwd(),'ksc.R'))               # k-Spectral Clustering (via Jungsuk)
source(file.path(getwd(),'timer.R'))             # adds tic() and toc() functions

hists = function(df,metric='sigma',zip='unspecified',norm=c()){
  .e <- environment() # capture local environment for use in ggplot
  dfsub = subset(df,select=c(metric,'id','model.name'))
  colnames(dfsub)[1] <- c('value')
  #dfm = melt(dfsub,id.vars=c('id'),variable.name='model.name')
  # plot several density plots at once color coded
  g = ggplot(dfsub, aes(x=value, color=model.name), environment=.e) + geom_density(size=1, alpha=0.2) + 
    ggtitle(paste(metric,'histograms for',zip))
  #png(paste("hist.png",sep=''),height=600,width=800)
  return(g)
}

bestFit = function(df,metric='rmse',sort=1,zip='unspecified') {
  inv = -1
  if(metric %in% c('rmse')) inv = 1
  n = dim(df[[metric]])[2]
  rnk = t(apply(inv*df[[metric]][,-1],MARGIN=1,FUN=rank))
  rnk = rnk[order(rnk[,sort]),]
  rnkm = melt(rnk)
  p = c()
  if (TRUE) {
    p = ggplot(rnkm, aes(x=X1,y=X2,fill=value)) +
      geom_tile() +
      scale_fill_gradient2("Rank",low='blue',mid="blue",high="grey",midpoint = 0) + # no key
      ylab("Model") +
      xlab("Residence") +
      coord_flip() +
      ggtitle(paste(metric,'fit ranking for',zip))
  } else {
    p = ggplot(rnkm, aes(value, color = X2)) + geom_bar(alpha=0.2) + 
       ggtitle(paste(metric,' rank histograms for',zip))
  }
  return(p) 
}

fitScatter = function(df,metric='rmse',zip='unspecified') {
  inv = -1
  if(metric %in% c('rmse')) inv = 1
  n = dim(df[[metric]])[2]
  rnk = t(apply(inv*df[[metric]][,-1],MARGIN=1,FUN=rank))
  rnk = rnk[order(rnk[,sort]),]
  rnkm = melt(rnk)
}

rankModel = function(m) {
  return(t(apply(m[,-1],MARGIN=1,FUN=rank)))
}

bestFitCount = function(df,metric='rmse',zip='unspecified') {
  inv = -1
  if(metric %in% c('rmse')) inv = 1
  n = dim(df[[metric]])[2]
  rnk = t(apply(inv*df[[metric]][,-1],MARGIN=1,FUN=rank))
  counts = apply(rnk<3,MARGIN=2,FUN=sum)
  barplot(counts,main=paste(metric,'count of 1st or 2nd rank for',zip),
          xlab="Model")
}


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
  
scalars = function(a,model.name=NULL,subset.name=NULL) {
  keepers = TRUE
  if(!is.null(model.name)) keepers = keepers & a$model.name==model.name
  if(!is.null(subset.name)) keepers = keepers & a$subset.name==subset.name
  b = subset(a,keepers,select=sapply(a,function(x) length(x[[1]]))==1) # select all columns with entries of length 1 (aka scalars)
  b$id = as.numeric(b$id) # if this is not numeric, the whole rbinded matrix is character and the data frame has factors instead of numeric cols
  # for some reason, all the data is stored as single entry list per data.frame cell this call pulls out
  # the actual valaues in those lists, so that the data.frame columns are the numeric arrays and factors we'd expect
  b$idx = 1:length(b$id) # add indices so we can do xxply operations on each row
  # for some reason, all the data is stored as single entry list per data.frame cell this call pulls out
  # the actual valaues in those lists, so that the data.frame columns are the numeric arrays and factors we'd expect
  c = ddply(b,.(idx),function(v) data.frame(lapply(v,function(x) x[[1]])))
  c$idx = c()
  return(c) # when above and this were one line, the cols were factors ?!!
}

cf = function(a,model.name,subset.name,col='Estimate') {
  keepers = (a$model.name==model.name & a$subset.name==subset.name)
  b = subset(a,keepers)
  allCols = data.frame(id=NA,t(b[[1,'contribution']])) # df with cols for all model params (including those that were aliased)
  b$id = as.numeric(b$id)                              # if this is not numeric, the whole rbinded matrix is character and the data frame has factors instead of numeric cols
  b$idx = 1:length(b$id)                               # setup ply to work on one row at a time
  c = dlply(b,.(idx),function(x) c(id=x$id,x$coefficients[[1]][,col])) # retrieve the named column from the coefficients matrix
  c$idx = c() # drop the index column
  g = as.numeric(lapply(c,length)) # calculate coefficient list lengths...
  c[g != Mode(g)] = c()            # remove the abnormal lengths caused by missing data
  h = as.numeric(lapply(c,function(x) any(x[-1] != 0))) # check for blank rows
  c[!h] = c() #remove rows of all zeros
  mergedDF = rbind.fill(allCols,data.frame(do.call(rbind,c))) # add missing cols using rbind.fill
  return(mergedDF[-1,]) # return the coefficient data, without the first extra row added to get the extra cols from rbind.fill
}
  
clean = function(a) {
  return (subset(a,a$total!=0))
}

demean = function(data,col_exp) {
  cols = grep(col_exp,colnames(data),value=TRUE)
  #row_var = apply(data[,cols],MARGIN=1,FUN=sd,na.rm=TRUE)
  row_mean = rowMeans(data[,cols],na.rm=TRUE)
  data[,cols] = (data[,cols] / row_mean)
  return(data)
}

p_summary = function(pvals,pmax=0.05,pattern='^HODWKWK') {
  n = dim(pvals)[1]
  m = dim(pvals)[2]
  pn = length(pattern)
  par(mfrow=c(pn,2))
  for (i in 1:pn) {
    significant = pvals[,-1] < pmax # boolean array of which values are stat. significant
    cols = grep(pattern[i],colnames(significant),value=TRUE)
    plot(colSums(significant[,cols],na.rm=TRUE), type='b',
         ylim=c(0,n),
         main=paste(pattern[i],' count p<',pmax,' per coefficient',sep=''),
         xlab='Coefficient',ylab='Count'  )
    hist(rowSums(significant[,cols],na.rm=TRUE),breaks=length(cols),
         #ylim=c(0,n),
         xlim=c(0,length(cols)+1),
         main=paste('Count hours p<',pmax,' per residence',sep=''),
         xlab='Hours'  )
  }
  par(mfrow=c(1,1))
}

c_summary = function(cvals,pvals,pmax=0.05,pattern='^HODWKWK') {
  significant = cbind(id=TRUE,pvals[,-1] < pmax) # boolean array of which values are stat. significant
  cvals[!significant] = NA
  n = dim(pvals)[1]
  m = dim(pvals)[2]
  pn = length(pattern)
  par(mfrow=c(pn,1))
  for (i in 1:pn) {
    cols = grep(pattern[i],colnames(significant),value=TRUE)
    plot(colMeans(cvals[,cols],na.rm=TRUE), type='b',
         #ylim=c(0,n),
         main=paste(pattern[i],' mean coef values p<',pmax,' N=',n,sep=''),
         xlab='Coefficient',ylab='coef / 24hr mean'  )
  }
  par(mfrow=c(1,1))
}


quad_chart = function(d,title) {
  
  # todo: update tehse for subsequent data runs:
  # the WK and ND suffixes have been moved to the beginning of the names
  # so "^HOD.+WK$" will become simply '^WKHOD'
  hodWK  = d[, c('id',grep("^HODWKWK", colnames(d), value=TRUE))]
  toutWK = d[, c('id',grep("^tout65.HODWKWK", colnames(d), value=TRUE))]
  hodND  = d[, c('id',grep("^HODWKND", colnames(d), value=TRUE))]
  toutND = d[, c('id',grep("^tout65.HODWKND", colnames(d), value=TRUE))]
  
  hodmWK =  melt(hodWK,id.vars=c('id'))
  toutmWK = melt(toutWK,id.vars=c('id'))
  hodmND =  melt(hodND,id.vars=c('id'))
  toutmND = melt(toutND,id.vars=c('id'))
  
  g1 = ggplot(hodmWK, aes(variable, value, group = id)) + geom_line(alpha = 0.05) + ylim(-0.1,3) + #ylim(-0.5,3) + 
    scale_x_discrete(labels=c('12am',1:11,'12pm',1:11)) + labs(title='Hour of day fixed effects (weekday)', x='Hour of day', y='fixed effect (kW)')
  g2 = ggplot(toutmWK, aes(variable, value, group = id)) + geom_line(alpha = 0.05) + ylim(-.01,4) + 
    scale_x_discrete(labels=c('12am',1:11,'12pm',1:11)) + labs(title='Hour of day temperature effects (weekday)', x='Hour of day', y='temperature effect (kW/deg above 65F)')
  g3 = ggplot(hodmND, aes(variable, value, group = id)) + geom_line(alpha = 0.05) + ylim(-0.1,3) + #ylim(-0.5,3) + 
    scale_x_discrete(labels=c('12am',1:11,'12pm',1:11)) + labs(title='Hour of day fixed effects (weekend)', x='Hour of day', y='fixed effect (kW)')
  g4 = ggplot(toutmND, aes(variable, value, group = id)) + geom_line(alpha = 0.05) + ylim(-.01,4) + 
    scale_x_discrete(labels=c('12am',1:11,'12pm',1:11)) + labs(title='Hour of day temperature effects (weekend)', x='Hour of day', y='temperature effect (kW/deg above 65F)')
  
  grid.arrange(g1,g2,g3,g4,nrow=2, as.table=FALSE, main=title)
}

combineSummaries = function(ziplist) {
  summaries = c()
  for (zip in ziplist) { 
    print(paste('loading data for',zip))
    load(file.path(getwd(),resultsDir,paste(zip,'_modelResults.RData',sep='')))
    summaries = rbind.fill(summaries,clean(modelResults$summaries))
  }
  return(summaries)
}

addZipData = function(summaries) {
  zipData = db.getZipData()
  summaries$zip = as.numeric(summaries$zip)
  return(merge(summaries,zipData,by.x='zip',by.y='zip5'))
}

model  = 'toutTOD_WKND'

allZips = c(94923,94503,94574,94559,94028,94539,94564,94702,94704,94085,
            95035,94041,95112,95113,95765,95648,95901,94531,94585,95205,
            95202,93619,93614,93304,93701,95631,95726,95223,95666)
summary = combineSummaries(allZips)

zip=93304
load(file.path(getwd(),resultsDir,paste(zip,'_modelResults.RData',sep='')))
summary   = clean(modelResults$summaries)
estimates = cf(summary,'toutPIECES24','all')
pvals     = cf(summary,'toutPIECES24','all',col='Pr(>|t|)')

slopeCols = grep('tout75_Inf',colnames(estimates),value=T)
hodCols   = grep('^HODH[0-9]+$',colnames(estimates),value=T)

slope75   = estimates[,c('id',slopeCols)]
slope75p  = pvals[,c('id',slopeCols)]
sig       = slope75p[,-1] > 0.05
slope75[cbind(rep(FALSE,dim(sig)[1]),sig)] <- NA
slope75m = melt(slope75,id.vars='id')
ggplot(data=slope75m,aes(x=variable, y=value,color=id)) + scale_x_discrete(labels=c('12am',1:11,'12pm',1:11)) + geom_line(aes(group=id))

hodConst  = estimates[,c('id',hodCols)]
hodConstp = pvals[,c('id',hodCols)]
sig       = hodConstp[,-1] > 0.05
#hodConst[cbind(rep(FALSE,dim(sig)[1]),sig)] <- NA
hodConstm = melt(hodConst,id.vars='id')
ggplot(data=hodConstm,aes(x=variable, y=value)) + scale_x_discrete(labels=c('12am',1:11,'12pm',1:11)) + geom_line(aes(group=id))

# create estimates for temperatures 0:100 by doing some linear algbra
pmCols    = grep('HODH15',colnames(estimates),value=T)
pmSlopes  = estimates[,c('id',pmCols)]
#pmConst   = estimates[,c('id',pmCols[1])]
pmSlopesp = pvals[,c('id',pmCols)]
X = regressor.piecewise(0:100,c(55,65,75)) # piecewise tout values from 0:100
X = cbind(1,X)                             # constant term
# one col of values per id
est = apply(as.matrix(pmSlopes[,-1]),MARGIN=1,FUN=function(row) t(X %*% as.matrix(row)))
colnames(est) <- pmSlopes[,1]
estm = melt(est)
ggplot(data=estm,aes(x=X1, y=value)) + geom_line(aes(group=X2))

#b = t(as.matrix(pmSlopes[5,-1]))
#plot(X %*% b)

zip = 94610
load(file.path(getwd(),resultsDir,paste(zip,'_modelResults.RData',sep='')))
summary   = clean(modelResults$summaries)
estimates.s = cf(summary,'toutTOD_WKND','summer')
estimates.w = cf(summary,'toutTOD_WKND','winter')
pvals.w     = cf(summary,'toutTOD_WKND','winter',col='Pr(>|t|)')
pvals.s     = cf(summary,'toutTOD_WKND','summer',col='Pr(>|t|)')

quad_chart(estimates.s,title=paste(zip,model,' (coast/moderate income', 'summer', 'only)'))


plot(colMeans(pvals.s[,grep("^HODWKWK", colnames(pvals.s), value=TRUE)]),main=paste(zip,' p-values'))

g1 = hists(scalars(summary,subset.name='summer'),'sigma',zip)
  
zip = 95223
load(file.path(getwd(),resultsDir,paste(zip,'_modelResults.RData',sep='')))
summary   = clean(modelResults$summaries)
estimates.s = cf(summary,model,'summer')
estimates.w = cf(summary,model,'winter')
pvals.w     = cf(summary,model,'winter',col='Pr(>|t|)')
pvals.s     = cf(summary,model,'summer',col='Pr(>|t|)')

# eliminate values not different from zero
pmax = 0.4 # maximum acceptable p-value in (0,1), with 0.05 being 95% confidence
estimates.s[cbind(FALSE,pvals.s[,-1] > pmax)] = NA
estimates.w[cbind(FALSE,pvals.w[,-1] > pmax)] = NA

# demean each set of coefficients
estimates.s = demean(estimates.s,"^tout65.HODWKWK")
estimates.w = demean(estimates.w,"^tout65.HODWKWK")
estimates.s = demean(estimates.s,"^tout65.HODWKND")
estimates.w = demean(estimates.w,"^tout65.HODWKND")
estimates.s = demean(estimates.s,"^HODWKWK")
estimates.w = demean(estimates.w,"^HODWKWK")
estimates.s = demean(estimates.s,"^HODWKND")
estimates.w = demean(estimates.w,"^HODWKND")

c_summary(estimates.s,pvals.s,0.05,pattern=c('^HODWKWK','^tout65.HODWKWK'))

p_summary(pvals.s,0.05,pattern=c('^HODWKWK','^tout65.HODWKWK'))


d_summary = clean(modelResults$d_summaries)

g1 = hists(scalars(d_summary,subset.name='all'),'r.squared',zip)
g2 = hists(scalars(d_summary,subset.name='all'),'sigma',zip)
grid.arrange(g1,g2)

# eliminate values not different from zero
#pmax = 0.05 # maximum acceptable p-value in (0,1), with 0.05 being 95% confidence
#estimates.s[cbind(FALSE,pvals.s[,-1] > pmax)] = NA
#estimates.w[cbind(FALSE,pvals.w[,-1] > pmax)] = NA

quad_chart(estimates.w,title=paste(zip,model,'(inland/lower income', 'winter', 'only)'))


par(mfrow=c(1,2))
plot(colMeans(pvals.w[,grep("^HODWKWK", colnames(pvals), value=TRUE)]),
     col='blue',
     ylim=c(0,1),
     main=paste(zip,'HOD p-values'),
     xlab='Hour of day',
     ylab='p-value (significant at < 0.05)')
points(colMeans(pvals.s[,grep("^HODWKWK", colnames(pvals), value=TRUE)]) / 2,pch=2,col='red')

plot(colMeans(pvals.w[,grep("^tout65.HODWKWK", colnames(pvals), value=TRUE)]),
     col='blue',
     ylim=c(0,1),
     main=paste(zip,'tout65xHOD p-values'),
     xlab='Hour of day',
     ylab='p-value (significant at < 0.05)')
points(colMeans(pvals.s[,grep("^tout65.HODWKWK", colnames(pvals), value=TRUE)]) / 2,pch=2,col='red')


g2 = hists(scalars(summary,subset.name='summer'),'sigma',zip)
grid.arrange(g2)
  
g2 = hists(scalars(summary,subset.name='all'),'r.squared',zip) + xlab('kWh/day')
grid.arrange(g2)


b = subset(modelResults$summaries,model.name == 'toutTOD' & subset.name == 'all',select=c(id,coefficients))
c = adply(b,.margins=1,.fun = function(row) {row$coefficients[[1]][1:10,'Estimate']})
c$coefficients = c()
c = c[, c('id',grep("^tout", colnames(c), value=TRUE))]
q = quantile(c$tout70_80,c(0.9))
c = subset(c,tout70_80 > q)
d = melt(c,id.vars=c('id'))
d$id = as.numeric(d$id) # somehow d$ids is a list of single entry lists?!
ggplot(d, aes(variable, value, group = id)) + geom_path(alpha = 0.1) + ylim(-0.25,0.5)

b = subset(modelResults$summaries,model.name == 'toutTOD' & subset.name == 'all',select=c(id,coefficients))
cnames = grep('^HOD',rownames(b[[1,'coefficients']]),value=TRUE)
c = adply(b,.margins=1,.fun = function(row) {row$coefficients[[1]][cnames,'Estimate']})
c$coefficients = c()
d = melt(c,id.vars=c('id'))
d$id = as.numeric(d$id) # somehow d$ids is a list of single entry lists?!
ggplot(d, aes(variable, value, group = id)) + geom_path(alpha = 0.05) #+ ylim(-0.1,0.15)

b = subset(modelResults$summaries,subset.name == 'all',select=c(id,model.name,rmse))
g = ggplot(b, aes(sigma, color = model.name)) + geom_density(size=1, alpha=0.2) + xlim(0,2.1)

a = subset(modelResults$features.basic,select=c(id,mean,min))
a$kWh = a$mean * 24 * 365
a$base = a$min * 24 * 365
b = subset(modelResults$summaries,model.name == 'toutPIECES' & subset.name == 'all',select=c(id,coefficients))
c = adply(b,.margins=1,.fun = function(row) {row$coefficients[[1]][,'Estimate']})
c$coefficients = c()
c = c[, c('id',grep("^tout", colnames(c), value=TRUE))]
q = quantile(c$tout70_80,c(0.75))
c$lower <- c$tout70_80 < q
c$quartile = paste('q',ceiling(rank(c$tout70_80) / (max(rank(c$tout70_80))/4)),seq='')
d = merge(a,c,by.y=c('id')) # get the kWh for the high Tout slope residence
d$rank = rank(d$kWh)
g = ggplot(a, aes(kWh)) + geom_density(size=1, alpha=0.2)
g = ggplot(d, aes(kWh,color=lower,fill=lower)) + geom_density(size=1, alpha=0.2) + geom_density(aes(kWh), size=1, alpha=0.2, color='black')
g = ggplot(d, aes(kWh,color=lower,fill=lower)) + geom_histogram(binwidth=500, alpha=.5, position="stack")
g = ggplot(d, aes(kWh,color=quartile,fill=quartile)) + geom_histogram(binwidth=500, alpha=.5, position="stack")
g = ggplot(d, aes(y=kWh,x=rank)) + geom_line()

resc[['z_94610']]$hourly$rmse[c('id','toutPIECES')]
gt = resc[['z_94610']]$hourly$rmse[c('id','toutPIECES')] < resc[['z_94610']]$hourly$rmse[c('id','standard')]
improved = gt[,2]     

                                                                                                    

