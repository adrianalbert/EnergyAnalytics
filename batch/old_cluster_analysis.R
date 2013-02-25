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
resultsDir = 'results' # real
#resultsDir = 'test_results' # sub -sample for testing

# run 'source' on all includes to load them 
source(file.path(conf.basePath,'localConf.R'))         # Sam's local computer specific configuration
source(file.path(conf.basePath,'ksc.R'))               # k-Spectral Clustering (via Jungsuk)
source(file.path(conf.basePath,'timer.R'))             # adds tic() and toc() functions

clusterPlots = function(regData,kscClusters){
  plot.list = list()
  for (i in 1:8){
    ggdata = melt(regData[(kscClusters$mem == i),],value.name='coef')
    if(is.null(ggdata$Var2)) {
      ggdata$Var1 = 1
      ggdata$Var2 = row.names(ggdata)
    }
    names(ggdata)[names(ggdata)=="Var2"] <- "TOD"
    names(ggdata)[names(ggdata)=="Var1"] <- "bldg"
    p = ggplot(ggdata, 
               aes(x=TOD, y=bldg, fill=coef)) + 
                 geom_tile() + 
                 scale_fill_gradient2("TOD fixed\neffect",guide=FALSE) + # no key
                 scale_x_discrete(breaks = 1:24) # no ticks or labels - not sure why
    plot.list[[i]] = p
  }
  return(plot.list)
}


RMSEHists = function(modelResults,zip='unspecified'){
  # convert from data.fram of lists to data.frame of numerics
  # TODO: this is suprisingly clumsy is there a cleaner way?
  
  RMSEData = data.frame( #MOY,
                         DOW=as.numeric(modelResults$summary[,'DOW']),
                         DOW_HOD=as.numeric(modelResults$summary[,'DOW_HOD']),
                         standard=as.numeric(modelResults$summary[,'standard']),
                         #HOW,
                         toutTOD=as.numeric(modelResults$summary[,'toutTOD']),
                         toutPIECES=as.numeric(modelResults$summary[,'toutPIECES'])
                         #toutPIECES_WKND  
                         )
  # turn into a data frame where the model type is 
  # just a text value in a column next to the RMSE (called 'variable')
  # This is the form that ggplot prefers
  RMSEStats = melt(RMSEData,variable.name='model',value.name='RMSE')
  # plot several density plots at once color coded
  g = ggplot(RMSEStats, aes(RMSE, color = model)) + geom_density(size=1, alpha=0.2) + 
             xlim(0, 0.6) + ggtitle(paste('RMSE histograms for',zip))
  #png(paste("hist.png",sep=''),height=600,width=800)
  return(g)
}


RMSEDiffs = function(modelResults,zip='unspecified'){
  # convert from data.fram of lists to data.frame of numerics
  # TODO: this is suprisingly clumsy is there a cleaner way?
  
  MOY=as.numeric(modelResults$summary[,'MOY'])
  DOW=as.numeric(modelResults$summary[,'DOW'])
  DOW_HOD=as.numeric(modelResults$summary[,'DOW_HOD'])
  standard=as.numeric(modelResults$summary[,'standard'])
  #HOW=as.numeric(modelResults$summary[,'HOW'])
  toutTOD=as.numeric(modelResults$summary[,'toutTOD'])
  toutPIECES=as.numeric(modelResults$summary[,'toutPIECES'])
  #toutPIECES_WKND=as.numeric(modelResults$summary[,'toutPIECES_WKND'])  

  RMSEDiffData = data.frame( 
    MOY = sort(DOW - MOY),
    standard = sort(DOW -standard),
    DOW_HOD = sort(DOW -DOW_HOD),
    #HOW = sort(DOW -HOW),
    toutTOD = sort(DOW -toutTOD),
    toutPIECES = sort(DOW -toutPIECES)
    #toutPIECES_WKND = sort(DOW -toutPIECES_WKND)
  )
  # turn into a data frame where the model type is 
  # just a text value in a column next to the RMSE (called 'variable')
  # This is the form that ggplot prefers
  RMSEDiffStats = melt(RMSEDiffData,variable.name='model',value.name='RMSEDiff')
  # plot several density plots at once color coded
  g = ggplot(data=RMSEDiffStats, aes(x=rep(1:length(MOY),5),y=RMSEDiff, colour=model)) + geom_line()
  #png(paste("hist.png",sep=''),height=600,width=800)
  return(g)
}

  # plot the count of group membership
  #cluster_table = table(kscClusters$mem)
  #hist(kscClusters$mem,length(cluster_table))
  
  #  p[[i]] = qplot(1:23,x=res_TOD$center[i,] * 100, data=data.frame(res_TOD$center),  geom="line", xlab='TOD', ylab='% of max', main=paste('C1 n=',TOD_table[i]))
#   p1 = qplot(1:length(kscClusters$center[1, ]),kscClusters$center[1,]   * 100, geom="line", xlab='TOD', ylab='% of max',main=paste('C0 n=', cluster_table[1]))
#   p2 = qplot(1:length(kscClusters$center[1, ]),kscClusters$center[2,]   * 100, geom="line", xlab='TOD', ylab='% of max',main=paste('C1 n=', cluster_table[2]))
#   p3 = qplot(1:length(kscClusters$center[1, ]),kscClusters$center[3,]   * 100, geom="line", xlab='TOD', ylab='% of max',main=paste('C2 n=', cluster_table[3]))
#   p4 = qplot(1:length(kscClusters$center[1, ]),kscClusters$center[4,]   * 100, geom="line", xlab='TOD', ylab='% of max',main=paste('C3 n=', cluster_table[4]))
#   p5 = qplot(1:length(kscClusters$center[1, ]),kscClusters$center[5,]   * 100, geom="line", xlab='TOD', ylab='% of max',main=paste('C4 n=', cluster_table[5]))
#   p6 = qplot(1:length(kscClusters$center[1, ]),kscClusters$center[6,]   * 100, geom="line", xlab='TOD', ylab='% of max',main=paste('C5 n=', cluster_table[6]))
#   p7 = qplot(1:length(kscClusters$center[1, ]),kscClusters$center[7,]   * 100, geom="line", xlab='TOD', ylab='% of max',main=paste('C6 n=', cluster_table[7]))
#   p8 = qplot(1:length(kscClusters$center[1, ]),kscClusters$center[8,]   * 100, geom="line", xlab='TOD', ylab='% of max',main=paste('C7 n=', cluster_table[8]))
#   grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8) #,p9,p10,p11,p12)


allResults = list.files(paste(conf.basePath,'/',resultsDir,sep=''),pattern = "[0-9]+_.+RData")
print(allResults)
for (resultFile in allResults) {
  zip = strsplit(resultFile,'_')[[1]][1]
  print(zip)
  if (zip != '93304') { next }
  
  # init the variables (e and modelResults) that load will assign values to
  e = c()               # any errors will loaded as e
  modelResults = list() # model run outcomes load as modelResults
  load(file.path(conf.basePath,resultsDir,paste(zip,'_modelResults.RData',sep='')))
  # TODO: we want to make a hist out of all the model results across every zipcode
  if(length(e)==0) { # the zip ran without a terminal error
    if(length(modelResults$ids)>0) {
      regData = modelResults$coef.HOW[,3:169]        # HOW
      regData = modelResults$coef.toutTOD[,26:192]   # HOW controlling toutTOD
      regData = modelResults$coef.toutPIECES[,-1:-8] # HOW controlling toutPIECES
      regData = modelResults$coef.standard[,-1:-7]   # TOD
      regData = modelResults$coef.toutPIECES[,2:8]   # toutPIECES
      regData = modelResults$coef.toutTOD[,2:25]     # toutTOD
      
      # normalize the results to get comparable values
      #regData = regData / modelResults$features.basic[,'97%']
      regData = regData[rowSums(is.na(regData))==0, ]
      kscClusters = ksc(regData,8)
      
      #print(modelResults)
      plot.list = list()
      plot.list[[1]] = RMSEHists(modelResults,zip)
      plot.list[[2]] = RMSEDiffs(modelResults,zip)
      plot.list[[3]] = clusterHist(kscClusters)
      plot.list = c(plot.list,clusterPlots(regData,kscClusters))
      do.call('grid.arrange', c(plot.list,list(ncol=2)))
      
      #RMSEHists(modelResults)
      break
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
