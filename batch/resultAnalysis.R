require(ggplot2)
require(gtools)
require(gridExtra)
require(reshape2)
require(plyr)
require(scales) # for muted

resultsDir = 'results_CP_solar' # real
#resultsDir = 'test_results' # sub -sample for testing

# run 'source' on all includes to load them 
source(file.path(getwd(),'localConf.R'))         # local computer specific configuration
source(file.path(getwd(),'ksc.R'))               # k-Spectral Clustering (via Jungsuk)
source(file.path(getwd(),'timer.R'))             # adds tic() and toc() functions

delist = function(df) {
  df$idx = 1:dim(df)[1]
  df2 <- ddply(df,'idx',function(X) apply(X,2,function(y) y[[1]])) # pull every column of every row of data out of its list of length 1
  df2 <- subset(df2, select = -c(idx) )
  return(df2)
}

scalars = function(a,model.name=NULL,subset.name=NULL) {
  keepers = TRUE
  if(!is.null(model.name)) keepers = keepers & a$model.name==model.name
  if(!is.null(subset.name)) keepers = keepers & a$subset.name==subset.name
  b = subset(a,keepers,select=sapply(a,function(x) (length(x[[1]]))==1 & class(x[[1]]) %in% c('numeric','character'))) # select all columns with entries of length 1 (aka scalars)
  b$id = as.numeric(b$id) # if this is not numeric, the whole rbinded matrix is character and the data frame has factors instead of numeric cols
  # for some reason, all the data is stored as single entry list per data.frame cell this call pulls out
  # the actual valaues in those lists, so that the data.frame columns are the numeric arrays and factors we'd expect
  return(delist(b)) # when above and this were one line, the cols were factors ?!!
}

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

hists = function(df,metric='sigma',zip='unspecified',norm=c()){
  .e <- environment() # capture local environment for use in ggplot
  dfsub = delist(subset(df,model.name != 'DOW', select=c(metric,'id','model.name')))
  colnames(dfsub)[1] <- c('value')
  dfsub$value = as.numeric(dfsub$value)
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



stderrs = function(a,model.name,subset.name=NA) {
  return(cf(a,model.name,subset.name,col='Std. Error'))
}

tvals = function(a,model.name,subset.name=NA) {
  return(cf(a,model.name,subset.name,col='t value'))
}

pvals = function(a,model.name,subset.name=NA) {
  return(cf(a,model.name,subset.name,col='Pr(>|t|)'))
}

cf = function(a,model.name,subset.name=NA,col='Estimate') {
  if (is.na(subset.name)) { 
    keepers = (a$model.name==model.name)
  }
  else {
    keepers = (a$model.name==model.name & a$subset.name==subset.name)
  }
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

combineSummaries = function(ziplist,resultType='summaries') {
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