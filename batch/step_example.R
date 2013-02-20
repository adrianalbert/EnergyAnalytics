# todo: is there a better way to detect the current directory?
conf.basePath = file.path('~/EnergyAnalytics/batch')
if(Sys.info()['sysname'] == 'Windows') {
  conf.basePath = file.path('f:/dev/pge_collab/EnergyAnalytics/batch')
} else {
  .libPaths('~/R/library') # use my local R library even from the comand line
}
setwd(conf.basePath)

# run 'source' on all includes to load them 
source(file.path(getwd(),'localConf.R'))         # Sam's local computer specific configuration
source(file.path(getwd(),'dbUtil.R'))            # database utilities
source(file.path(getwd(),'DataClasses.R'))       # Object code for getting meter and weather data 
source(file.path(getwd(),'ksc.R'))               # k-Spectral Clustering (via Jungsuk)
source(file.path(getwd(),'basicFeatures.R'))     # typical max, min, mean, range
source(file.path(getwd(),'regressionSupport.R')) # mostly regressor manipulation
source(file.path(getwd(),'timer.R'))             # adds tic() and toc() functions

r = ResDataClass(820735863,94610)
dfd = regressorDFAggregated(r,rm.na=TRUE)$daily
df  = regressorDF(r,rm.na=TRUE)
full = "kwh ~ tout.mean + pout.mean + rh.mean + WKND + vac"
toutTOD = 'kw ~ 0 + tout65:HODWK + HODWK'
m0   = lm(kwh ~  1,dfd)
mfull = lm(kwh ~ tout.mean + pout.mean + rh.mean + WKND + vac,dfd)
#stepped = step(m0,direction='forward',k=1,
#                  scope=kwh ~ tout.mean + pout.mean + rh.mean + WKND + vac)
#as = stepped$anova
#ss = summary(stepped)



stepped2 = step(lm(kw ~ 0,df),direction='forward',k=2,
                scope=kw ~ 0 + tout65:HOD + pout + rh + tout65_d1 + HOD)
#               scope=kw ~ 0 + tout65:HODWK + HODWK)
as2 = stepped2$anova
ss2 = summary(stepped2)


lm(formula = kw ~ HOD + rh + pout + tout65_d1 + HOD:tout65 - 
     1, data = df)

