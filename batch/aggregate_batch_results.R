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

conf.dataPath = 'F:/dev/pge_collab/EnergyAnalytics/batch'

# resultsDir = 'results' # real
resultsDir = 'results' # sub -sample for testing

# run 'source' on all includes to load them 
source(file.path(conf.basePath,'localConf.R'))         # Sam's local computer specific configuration
source(file.path(conf.basePath,'ksc.R'))               # k-Spectral Clustering (via Jungsuk)
source(file.path(conf.basePath,'timer.R'))             # adds tic() and toc() functions

allResults = list.files(paste(conf.dataPath,resultsDir,sep=''),pattern = "[0-9]+_.+RData")
allSummary = data.frame()
idx = 0
n = length(allResults)
for (resultFile in allResults) {
  idx = idx + 1
  zip = strsplit(resultFile,'_')[[1]][1]
  print(paste(zip,' (',idx,'/',n,')'))
  #if (zip != '94610') { next }
  
  # init the variables (e and modelResults) that load will assign values to
  e = c()               # any errors will loaded as e
  modelResults = list() # model run outcomes load as modelResults
  load(file.path(conf.dataPath,resultsDir,paste(zip,'_modelResults.RData',sep='')))
  # TODO: we want to make a hist out of all the model results across every zipcode
  if(length(e)==0) { # the zip ran without a terminal error
    if(length(modelResults$ids)>0) {
      allSummary = rbind(allSummary, modelResults$features.basic[,c('mean','min')])
      
      #print(modelResults)
    } else {
      print(paste('No results for',zip))
    }
  }
  if(length(e)>1) { # the zip had a problem if e is a list...
    print(paste(zip,'had a problem:',e$message,'on call to:'))
    print(e$call)
  }
  rm(e)
  #rm(modelResults) 
}
allSummary = allSummary * 365*24
colnames(allSummary) <- c('kWh','kWh.base')

plot(cumsum(allSummary[with(allSummary, order(-kWh)), 'kWh'])/1000,xlab='index',ylab='MWh/yr')
plot(allSummary[,'kWh.base'],allSummary[,'kWh'])
hist(allSummary[,'kWh'],breaks=500,xlim=c(0,40000),ylab='count',xlab='Annual kWh',main=paste('Annual energy use (kWh/yr) for n =',dim(allSummary)[1]))

