#.libPaths('~/R/library') # use my local R library even from the command line

require(ggplot2)
require(gtools)
require(gridExtra)
require(reshape2)
require(plyr)
require(scales) # for muted

# todo: is there a better way to detect the current directory?
conf.basePath = file.path('~/EnergyAnalytics/batch')
if(Sys.info()['sysname'] == 'Windows' & Sys.info()['user'] == 'Sam') {
  conf.basePath = file.path('f:/dev/pge_collab/EnergyAnalytics/batch')
}
setwd(conf.basePath)
source(file.path(getwd(),'resultAnalysis.R'))
source(file.path(getwd(),'zipMap.R'))

resultsDir = 'results_daily' # real

if(exists('cfg')) { resultsDir = cfg$outDir }

allZips = c(94923,94503,94574,94559,94028,94539,94564,94702,94704,94085,
            95035,94041,95112,95113,95765,95648,95901,94531,94585,95205,
            95202,93619,93614,93304,93701,95631,95726,95223,95666)
#summary = combineSummaries(allZips,resultType='d_summaries') # fails from too much data?!

allZips = c(94923,94503,94574,94559)
#comb = combineSummaries(allZips,resultType='d_summaries')
basics     = combine(allZips,resultType='features.basic',as.matrix,appendZipData=T)
basicMeans = combine(allZips,resultType='features.basic',function(x) { t(colMeans(x)) },appendZipData=T )
sclrs      = combine(allZips,resultType='d_summaries',scalars,appendZipData=T)

cfs  = list()
stde = list()
pvs  = list()
tvs  = list()
models = unique(sclrs$model.name)
for(model in models) {
  cfs[[model]]  = combine(allZips,resultType='d_summaries',cf,     model.name=model,appendZipData=T)
  stde[[model]] = combine(allZips,resultType='d_summaries',stderrs,model.name=model,appendZipData=T)
  pvs[[model]]  = combine(allZips,resultType='d_summaries',pvals,  model.name=model,appendZipData=T)
  tvs[[model]]  = combine(allZips,resultType='d_summaries',tvals,  model.name=model,appendZipData=T)
}

# find the best r.squared result for each sp_id
a = do.call(rbind,by(sclrs,sclrs$id,function(df) df[which.max(df$r.squared),]))

hists(sclrs,metric='r.squared')
sclrsZ = addZipData(sclrs)
zipHists(sclrsZ,metric='r.squared',model.name='toutDailyCP')


zip='94923'
#d_summary = modelResults$d_summary #combineSummaries(zips,'d_summaries')
load(file.path(getwd(),resultsDir,paste(zip,'_modelResults.RData',sep='')))
summary   = clean(modelResults$d_summaries)
cp = modelResults$changePoints
cpgrid = t(apply(cp,1,function(X) X[['changePoints']]['cp',]))
changes = ksc(cpgrid,8)
plotClusterCenters(changes)
showClusters(changes)
clusterHist(changes)
p = clusterHeatMaps(cpgrid,changes)
do.call(grid.arrange,p)

hists(summary,metric='sigma',zip=zips[2])
hists(summary,metric='r.squared',zip=zips[2])
hists(summary,metric='adj.r.squared',zip=zips[2])

hists(d_summary,metric='sigma',zip=zips[1])
hists(d_summary,metric='r.squared',zip=zips[1])
hists(d_summary,metric='adj.r.squared',zip=zips[1])

estimates = cf(summary,'toutPIECES24','all')
pvals     = cf(summary,'toutPIECES24','all',col='Pr(>|t|)')

slopeCols = grep('tout75_Inf',colnames(estimates),value=T)[-c(1:8)]
slope75   = estimates[,c('id',slopeCols)]
slope75p  = pvals[,c('id',slopeCols)]
sig       = slope75p[,-1] >= 0.0
slope75[cbind(rep(FALSE,dim(sig)[1]),!sig)] <- NA
slope75m = melt(slope75,id.vars='id')
ggplot(data=slope75m,aes(x=variable, y=value)) + scale_x_discrete(labels=c('12am',1:11,'12pm',1:11)) + geom_line(aes(group=id),alpha = 0.05)

a = as.matrix(slope75[,-1])



hodCols   = grep('^HODH[0-9]+$',colnames(estimates),value=T)
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
colnames(est) <- paste('id_',pmSlopes[,1],sep='') # restore id names
estm = melt(est)
ggplot(data=estm,aes(x=X1, y=value,color=X2)) + geom_line(aes(group=X2))


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

                                                                                                    

