#!/usr/bin/Rscript

#.libPaths('~/R/library') # use my local R library even from the comand line

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
resultsDir = 'results_tout_schedule' # real
#resultsDir = 'test_results' # sub -sample for testing

# run 'source' on all includes to load them 
source(file.path(conf.basePath,'localConf.R'))         # Sam's local computer specific configuration
source(file.path(conf.basePath,'ksc.R'))               # k-Spectral Clustering (via Jungsuk)
source(file.path(conf.basePath,'timer.R'))             # adds tic() and toc() functions

consolidateResults = function(modelResults) {
  
  # TODO: set up a data structures for these more complex outputs
  #modelResults$summaries[[1,'DOW']]$coefficients
  #modelResults$summaries[[1,'DOW']]$contribution
  
  res = list(
      features = modelResults$features.basic,
      hourly   = metricSummary(modelResults$summaries),
      daily    = metricSummary(modelResults$d_summaries),
      monthly  = metricSummary(modelResults$m_summaries)
    )
  return(res)
}

toRow = function(res,var) {
  return(data.frame(lapply(res,function(x) x[[var]])))
}

resultDF = function(modelGroup,var) {
  res = adply(modelGroup,.margins=1,.fun=toRow,var=var)
  res$X1 <- c() # some sort of artifact of adply
  return(res)
}

metricSummary = function(modelGroup) {
  return( list (
    rmse          = cbind(id=modelResults$id,resultDF(modelGroup,'sigma')),
    r.squared     = cbind(id=modelResults$id,resultDF(modelGroup,'r.squared')),
    adj.r.squared = cbind(id=modelResults$id,resultDF(modelGroup,'adj.r.squared')),
    log.liklihood = cbind(id=modelResults$id,resultDF(modelGroup,'logLik')),
    aic           = cbind(id=modelResults$id,resultDF(modelGroup,'AIC')),
    total         = cbind(id=modelResults$id,resultDF(modelGroup,'total'))
  ))
}

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

allResults = list.files(paste(conf.basePath,'/',resultsDir,sep=''),pattern = "[0-9]+_modelResults.RData")
resr = list()
resc = list()
for (resultFile in allResults) {
  zip = strsplit(resultFile,'_')[[1]][1]
  print(zip)
  if (zip %in% c(93727)) { next }
  e = c()
  load(file.path(conf.basePath,resultsDir,paste(zip,'_modelResults.RData',sep='')))
  if(length(e)==0) { # the zip ran without a terminal error
    resc[[paste('z_',zip,sep='')]] = consolidatedResults
    resr[[paste('z_',zip,sep='')]] = modelResults
    #consolidatedResults = consolidateResults(modelResults)
    #save(consolidatedResults, file=file.path(conf.basePath,resultsDir,paste(zip,'_consolidatedResults.RData',sep='')))
    print(names(consolidatedResults))
    print(names(modelResults))
    # "rmse"          "r.squared"     "adj.r.squared" "log.liklihood" "aic" 
    # make a heat map of the goodness of fit
    #g = bestFit(res$hourly,'rmse',sort='standard',zip=zip); g
    #p = hists(res$hourly,'rmse',zip=zip); p
    # make a line plot of the goodness of fit    
  }
  if(length(e)>1) { # the zip had a problem if e is a list...
    print(paste(zip,'had a problem:',e$message,'on call to:'))
    print(e$call)
  }
  rm(e)
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

cf = function(a,model.name,subset.name) {
  keepers = (a$model.name==model.name & a$subset.name==subset.name)
  b = subset(a,keepers)
  b$id = as.numeric(b$id) # if this is not numeric, the whole rbinded matrix is character and the data frame has factors instead of numeric cols
  b$idx = 1:length(b$id)
  c = dlply(b,.(idx),function(x) c(id=x$id,x$coefficients[[1]][,'Estimate']))
  c$idx = c()
  g = as.numeric(lapply(c,length)) # calculate coefficient list lengths...
  c[g != Mode(g)] = c()            # remove the abnormal lengths
  h = as.numeric(lapply(c,function(x) any(x[-1] != 0))) # check for blank rows
  c[!h] = c() #remove rows of all zeros
  return(data.frame(do.call(rbind,c))) # when above and this were one line, the cols were factors ?!!
}
  
clean = function(a) {
  return (subset(a,a$total!=0))
}



quad_chart = function(d,title) {
  hodWK  = d[, c('id',grep("^HOD.+WK$", colnames(d), value=TRUE))]
  toutWK = d[, c('id',grep("^tout65.+WK$", colnames(d), value=TRUE))]
  hodND  = d[, c('id',grep("^HOD.+ND$", colnames(d), value=TRUE))]
  toutND = d[, c('id',grep("^tout65.+ND$", colnames(d), value=TRUE))]
  
  hodmWK =  melt(hodWK,id.vars=c('id'))
  toutmWK = melt(toutWK,id.vars=c('id'))
  hodmND =  melt(hodND,id.vars=c('id'))
  toutmND = melt(toutND,id.vars=c('id'))
  
  g1 = ggplot(hodmWK, aes(variable, value, group = id)) + geom_line(alpha = 0.05) + ylim(-1,3) + 
    scale_x_discrete(labels=c('12am',1:11,'12pm',1:11)) + labs(title='Hour of day fixed effects (weekday)', x='Hour of day', y='fixed effect (kW)')
  g2 = ggplot(toutmWK, aes(variable, value, group = id)) + geom_line(alpha = 0.05) + ylim(-.1,0.2) + 
    scale_x_discrete(labels=c('12am',1:11,'12pm',1:11)) + labs(title='Hour of day temperature effects (weekday)', x='Hour of day', y='temperature effect (kW/deg above 65F)')
  g3 = ggplot(hodmND, aes(variable, value, group = id)) + geom_line(alpha = 0.05) + ylim(-1,3) + 
    scale_x_discrete(labels=c('12am',1:11,'12pm',1:11)) + labs(title='Hour of day fixed effects (weekend)', x='Hour of day', y='fixed effect (kW)')
  g4 = ggplot(toutmND, aes(variable, value, group = id)) + geom_line(alpha = 0.05) + ylim(-.1,0.2) + 
    scale_x_discrete(labels=c('12am',1:11,'12pm',1:11)) + labs(title='Hour of day temperature effects (weekend)', x='Hour of day', y='temperature effect (kW/deg above 65F)')
  
  grid.arrange(g1,g2,g3,g4,nrow=2, as.table=FALSE, main=title)
}

zip = 94610
load(file.path(conf.basePath,resultsDir,paste(zip,'_modelResults.RData',sep='')))
#d = cf(modelResults$summaries,'toutTOD','summer')
quad_chart(cf(clean(modelResults$summaries),'toutTOD_WKND','summer'),title=paste(zip,' (coast/moderate $) tout65:HOD + HOD model (summer only)'))

g1 = hists(scalars(clean(modelResults$summaries),subset.name='all'),'sigma',zip)
  
zip = 93304
load(file.path(conf.basePath,resultsDir,paste(zip,'_modelResults.RData',sep='')))
#d = cf(modelResults$summaries,'toutTOD','summer')
quad_chart(cf(clean(modelResults$summaries),'toutTOD_WKND','summer'),title=paste(zip,' (inland/moderate $) tout65:HOD + HOD model (summer only)'))

g2 = hists(scalars(clean(modelResults$summaries),subset.name='all'),'r.squared',zip)
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

                                                                                                    

