require(ggplot2)
require(gtools)
require(gridExtra)
require(reshape2)
require(plyr)
require(scales) # for muted

# run 'source' on all includes to load them 
source(file.path(getwd(),'localConf.R'))         # local computer specific configuration
source(file.path(getwd(),'ksc.R'))               # k-Spectral Clustering (via Jungsuk)
source(file.path(getwd(),'timer.R'))             # adds tic() and toc() functions


# Helper function ---------------------------------------------------------

delist = function(df) {
  df2 = df
  for(col in colnames(df)) {
    df2[col] = unlist(df[col])
  }
  #df$idx = 1:dim(df)[1]
  #df2 <- ddply(df,'idx',function(X) apply(X,2,function(y) y[[1]])) # pull every column of every row of data out of its list of length 1
  #df2 <- subset(df2, select = -c(idx) )
  return(df2)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
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

# Run data extraction -----------------------------------------------------

scalars = function(a,model.name=NULL,subset.name=NULL) {
  keepers = TRUE
  if(!is.null(model.name)) keepers = keepers & a$model.name==model.name
  if(!is.null(subset.name)) keepers = keepers & a$subset.name==subset.name
  b = subset(a,keepers,select=sapply(a,function(x) (length(x[[1]]))==1 & class(x[[1]]) %in% c('numeric','character'))) # select all columns with entries of length 1 (aka scalars)
  return(delist(b))
}

cf = function(a,model.name,subset.name=NULL,col='Estimate',cfName='coefficients') {
  sn = T
  mn = T
  #print(names(a))
  if(is.null(cfName)) { cfName='coefficients' }
  mn = unlist(a[,'model.name'])==model.name
  if (! is.null(subset.name)) { 
    sn = unlist(a[,'subset.name'])==subset.name
  }
  b = subset(a,mn & sn) # filter to the subset and model requested
  if(empty(b)) { return(c()) }
  out = t(apply(b,1,function(x) return(c(id=x$id,x[cfName][[1]][,col]))))
  if(class(out[1]) == 'list') { # this seems to happen when some regessons had fewer non-NA coefficients than others.
    print('warning. Some regressons have incomplete coefficients and are being dropped.')
    lengths = as.numeric(lapply(out,length))
    out = do.call(rbind,out[lengths == Mode(lengths)])
  }
  return(out)
}

pvalsNW    = function(a,model.name,subset.name=NULL) {
  return(pvals(a,model.name,subset.name,cfName='coefficientsNW'))
}

stderrsNW = function(a,model.name,subset.name=NULL) {
  return(stderrs(a,model.name,subset.name,cfName='coefficientsNW'))
}

stderrs = function(a,model.name,subset.name=NULL,cfName=NULL) {
  return(cf(a,model.name,subset.name,col='Std. Error',cfName=cfName))
}

tvals = function(a,model.name,subset.name=NULL,cfName=NULL) {
  return(cf(a,model.name,subset.name,col='t value',cfName=cfName))
}

pvals = function(a,model.name,subset.name=NULL,cfName=NULL) {
  return(cf(a,model.name,subset.name,col='Pr(>|t|)',cfName=cfName))
}


# Aggregation function ----------------------------------------------------

combineSummaries = function(ziplist,resultType='summaries') {
  summaries = c()
  i = 0
  n = length(ziplist)
  for (zip in ziplist) { 
    i = i+1
    print(paste('loading data for',zip,'(',i,'/',n,')'))
    dataFile = file.path(getwd(),resultsDir,paste(zip,'_modelResults.RData',sep=''))
    if (! file.exists(dataFile)){
      print(paste('No data file for',zip,'skipping.'))
      next
    }
    load(dataFile)
    summaries = rbind(summaries,clean(modelResults[[resultType]]))
    rm(modelResults)
  }
  return(summaries)
}

combine = function(ziplist,resultType='summaries',subResultType=NULL,fun=function(x) { x },model.name=NULL,subset.name=NULL,appendZipData=F) {
  result = c()
  rList = as.list(rep(NA,600)) # there are < 600 zips so far... 
  i = 0
  n = length(ziplist)
  prevColCount <- 0
  for (zip in ziplist) { 
    i = i+1
    zip = as.numeric(zip)
    print(paste('loading data for',zip,'(',i,'/',n,')'))
    dataFile = file.path(getwd(),resultsDir,paste(zip,'_modelResults.RData',sep=''))
    #print(dataFile)
    if (! file.exists(dataFile)){
      print(paste('No data file for',zip,'skipping.'))
      next
    }
    load(dataFile)
    if(! resultType %in% names(modelResults) ) {
      print(paste('no entry for',resultType))
      next
    }
    data = modelResults[[resultType]]
    if(! is.null(subResultType)) {
      data = data[[subResultType]]
      if(is.null(data)) {
        print(paste('WARNING. No data found under subResultType',subResultType))
      }
    }
    if(length(data) == 0) {
      print(paste('no results for',zip,'. Skipping.'))
      next
    }
    #print(class(modelResults[['d_summaries']]))
    #print(class(modelResults$d_summaries))
    #print(cf(modelResults[['d_summaries']],model.name=model.name,subset.name=NULL))
        
    if(is.null(model.name)) {
      new = fun(data)
    } else {
      new = fun(data,model.name=model.name,subset.name=subset.name)
    }
    
    if(length(new) == 0) { next }
    # the next several lines try to eliminate duplicate entries
    dataNames = names(new) # assume data.frame
    if(is.null(dataNames)) dataNames = colnames(new) # fall back to array with names cols
    #print(dataNames)
    uniqueCols = intersect(dataNames,c('id','model.name','subset.name')) # find all cols of interest
    if(length(uniqueCols) > 0) {
      new = new[!duplicated(new[,uniqueCols]),] # somehow, some of the sp_id entries are duplicates
    }
    # cbind doesn't work on a single row...
    if(is.null(dim(new))) {
      new = matrix(new,nrow=1,dimnames=list(NULL,names(new)))
    }
    new = cbind(new,zip5=zip) # ensure the zipcode is there
    newLength = length(colnames(new))
    if(prevColCount > newLength) {
      print(colnames(new))
      print('trouble with column length')
      print(new)
      # skip?
    }
    prevColCount <- newLength
    rownames(new) <- c()
    rList[[i]] = new # inserting into a pre-allocated list is faster than running rbind all the time
    #rList[[length(rList)+1]] = new
    rm(modelResults)
  }
  #print(rList[[38]])
  #print(colnames(rList[[38]]))
  nCols = lapply(rList,function(x) length(colnames(x)))
  rList[nCols == 0] <- NULL # remove empty parts of the list.
  #print(do.call(rbind,lapply(rList,function(x) c(length(colnames(x)),unique(x['zip5'])))))
  #print(do.call(rbind,lapply(rList,function(x) c(colnames(x),unique(x['zip5'])))))
  
  result = do.call(rbind,rList[!is.na(rList)]) # here we bind the results together
  if(appendZipData) { result = addZipData(result) }
  return(result)
}

addZipData = function(orig) {
  if(! "data.frame" %in% class(orig)) {  # because the zip data is mixed, the return value
    orig = data.frame(orig)              # has to be a data frame too
  }
  zipData = DATA_SOURCE$getZipData(useCache=T)
  zipCol = grep('^zip',colnames(orig),value=T)[1]
  print(colnames(orig))
  print('1')
  print(zipCol)
  orig[[zipCol]] = as.numeric(orig[[zipCol]])
  return(merge(orig,zipData,by.x=zipCol,by.y='zip5'))
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

# Plotting function -------------------------------------------------------

clusterHist = function(kscClusters) {
  df = as.data.frame(list(cluster=kscClusters$mem))
  maxX= length(unique(kscClusters$cluster))
  p = qplot(cluster, data=df, geom="histogram",binwidth=0.5) +
    scale_x_discrete(breaks=0:8) + ggtitle('Count of cluster members')
  return(p)
}

plotClusterCenters = function(kscClusters) {
  clusterCounts = table(kscClusters$mem)
  p = list()
  n = dim(kscClusters$center)[1]
  par(mfrow=c(2,ceiling(n/2)))
  for (i in 1:n) {
    values = kscClusters$center[i,]
    
    plot(values,type='b',xlab='TOD',ylab='change point (normalized)',
         main=paste('C',i,' n=',clusterCounts[i],sep=''))
    #p[[i]] = qplot(1:length(values),values, 
    #               geom='line',
    #               xlab='TOD',ylab='change point (F)',
    #               main=paste('C',i,' n=',clusterCounts[i]))
  }
  #print(p)
  #do.call(grid.arrange,p)
}

showClusters = function(kscClusters) {
  df = data.frame(kscClusters$center)
  df$idx = 1:dim(kscClusters$center)[1]
  dfm = melt(df,id.vars=c('idx'))
  ggplot(dfm,aes(x=variable,y=value)) + geom_line(aes(group=idx)) 
}

clusterHeatMaps = function(regData,kscClusters){
  clusterCounts = table(kscClusters$mem)
  plot.list = list()
  for (i in 1:max(kscClusters$mem)){
    ggdata = melt(regData[(kscClusters$mem == i),])
    if(is.null(ggdata$X2)) {
      ggdata$X1 = 1
      ggdata$X2 = row.names(ggdata)
    }
    names(ggdata)[names(ggdata)=="X2"] <- "hour"
    names(ggdata)[names(ggdata)=="X1"] <- "idx"
    p = ggplot(ggdata, 
               aes(x=hour, y=idx, fill=value)) + 
      geom_tile() + 
      scale_fill_gradient2("change point",midpoint=70,low="red",mid="white",high="blue") + # no key
      scale_x_discrete(breaks = 1:24) + # no ticks or labels - not sure why
      ggtitle(paste('C',i,' n=',clusterCounts[i],sep=''))
    plot.list[[i]] = p
  }
  return(plot.list)
}

hists = function(df,metric='sigma',zip='unspecified',norm=c(),xlim=NULL,ylim=NULL,xlab=NULL,ylab=NULL){
  .e <- environment() # capture local environment for use in ggplot
  dfsub = delist(subset(df,model.name != 'DOW', select=c(metric,'id','model.name')))
  colnames(dfsub)[1] <- c('value')
  dfsub$value = as.numeric(dfsub$value)
  #dfm = melt(dfsub,id.vars=c('id'),variable.name='model.name')
  # plot several density plots at once color coded
  g = ggplot(dfsub, aes(x=value, color=model.name), environment=.e) + geom_density(size=1, alpha=0.2) + 
    labs(title=paste(metric,'density for',zip))
  if(! is.null(xlim)) { g = g + xlim(xlim[1],xlim[2]) }
  if(! is.null(ylim)) { g = g + ylim(ylim[1],ylim[2]) }
  #png(paste("hist.png",sep=''),height=600,width=800)
  return(g)
}

zipHists = function(df,metric='sigma',model.name,norm=c(),zip=NULL,xlim=NULL,ylim=NULL,xlab=NULL,ylab=NULL){
  .e <- environment() # capture local environment for use in ggplot
  filter = df$model.name == model.name 
  if(! is.null(zip)) { filter = filter & df$zip5==zip }
  dfsub = delist(subset(df,filter,select=c(metric,'id','zip5','cecclmzn')))
  colnames(dfsub)[1] <- c('value')
  dfsub$value = as.numeric(dfsub$value)
  #dfm = melt(dfsub,id.vars=c('id'),variable.name='model.name')
  # plot several density plots at once color coded
  g = ggplot(dfsub, aes(x=value, color=cecclmzn), environment=.e) + geom_density(size=1, alpha=0.2) + 
    labs(title=paste(metric,'density for',model.name)) 
  if(! is.null(xlim)) { g = g + xlim(xlim[1],xlim[2]) }
  if(! is.null(ylim)) { g = g + ylim(ylim[1],ylim[2]) }
  if(! is.null(xlab)) { g = g + labs(x=xlab) }
  if(! is.null(ylab)) { g = g + labs(y=ylab) }
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
  hodWK  = d[, c('id',grep("^HODWKWK",        colnames(d), value=TRUE))]
  toutWK = d[, c('id',grep("^tout65.HODWKWK", colnames(d), value=TRUE))]
  hodND  = d[, c('id',grep("^HODWKND",        colnames(d), value=TRUE))]
  toutND = d[, c('id',grep("^tout65.HODWKND", colnames(d), value=TRUE))]
  
  hodmWK  = melt(hodWK, id.vars=c('id'))
  toutmWK = melt(toutWK,id.vars=c('id'))
  hodmND  = melt(hodND, id.vars=c('id'))
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

