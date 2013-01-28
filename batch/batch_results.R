#!/usr/bin/Rscript

#.libPaths('~/R/library') # use my local R library even from the comand line

require(ggplot2)
require(gtools)
require(gridExtra)
require(reshape2)

# todo: is there a better way to detect the current directory?
conf.basePath = file.path('~/EnergyAnalytics/batch')
if(Sys.info()['sysname'] == 'Windows') {
  conf.basePath = file.path('f:/dev/pge_collab/EnergyAnalytics/batch')
}
resultsDir = 'results_summer' # real
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

hists = function(df,metric='rmse',zip='unspecified',norm=c()){
  # This is the form that ggplot prefers
  raw = df[[metric]]
  if(length(norm)>0) {
    raw = raw / norm
  }
  dfm = melt(raw,id.vars=c('id'),variable.name='model')
  # plot several density plots at once color coded
  g = ggplot(dfm, aes(value, color = variable)) + geom_density(size=1, alpha=0.2) + 
    xlim(0, 1) + ggtitle(paste(metric,'histograms for',zip))
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
  #countsm = melt(counts)
  #p = ggplot(countsm, aes(value)) +
  #  geom_bar() +
  #  ylab("Count") +
  #  xlab("Model") +
  #  ggtitle(paste(metric,'fit ranking for',zip))
  #return(p) 
}

# allResults = list.files(paste(conf.basePath,'/',resultsDir,sep=''),pattern = "[0-9]+_.+RData")
# res = c()
# for (resultFile in allResults) {
#   zip = strsplit(resultFile,'_')[[1]][1]
#   print(zip)
#   if (zip %in% c(93727)) { next }
#   e = c()
#   res = c()
#   load(file.path(conf.basePath,resultsDir,paste(zip,'_modelResults.RData',sep='')))
#   if(length(e)==0) { # the zip ran without a terminal error
#     consolidatedResults = consolidateResults(modelResults)
#     save(consolidatedResults, file=file.path(conf.basePath,resultsDir,paste(zip,'_consolidatedResults.RData',sep='')))
#     print(names(consolidatedResults))
#     # "rmse"          "r.squared"     "adj.r.squared" "log.liklihood" "aic" 
#     # make a heat map of the goodness of fit
#     #g = bestFit(res$hourly,'rmse',sort='standard',zip=zip); g
#     #p = hists(res$hourly,'rmse',zip=zip); p
#     # make a line plot of the goodness of fit    
#   }
#   if(length(e)>1) { # the zip had a problem if e is a list...
#     print(paste(zip,'had a problem:',e$message,'on call to:'))
#     print(e$call)
#   }
#   rm(e)
# }

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

b = subset(modelResults$summaries,model.name == 'toutPIECES' & subset.name == 'all',select=c(id,coefficients))
c = adply(b,.margins=1,.fun = function(row) {row$coefficients[[1]][1:10,'Estimate']})
c$coefficients = c()
c = c[, c('id',grep("^tout", colnames(c), value=TRUE))]
q = quantile(c$tout70_80,c(0.9))
c = subset(c,tout70_80 > q)
d = melt(c,id.vars=c('id'))
d$id = as.numeric(d$id) # somehow d$ids is a list of single entry lists?!
ggplot(d, aes(variable, value, group = id)) + geom_path(alpha = 0.1) + ylim(-0.25,0.5)

b = subset(modelResults$summaries,model.name == 'toutTOD' & subset.name == 'all',select=c(id,coefficients))
cnames = grep('^tout',rownames(b[[1,'coefficients']]),value=TRUE)
c = adply(b,.margins=1,.fun = function(row) {row$coefficients[[1]][cnames,'Estimate']})
c$coefficients = c()
d = melt(c,id.vars=c('id'))
d$id = as.numeric(d$id) # somehow d$ids is a list of single entry lists?!
ggplot(d, aes(variable, value, group = id)) + geom_path(alpha = 0.2) + ylim(-0.1,0.15)

b = subset(modelResults$summaries,subset.name == 'all',select=c(id,model.name,rmse))



resc[['z_94610']]$hourly$rmse[c('id','toutPIECES')]
gt = resc[['z_94610']]$hourly$rmse[c('id','toutPIECES')] < resc[['z_94610']]$hourly$rmse[c('id','standard')]
improved = gt[,2]     

                                                                                                    

