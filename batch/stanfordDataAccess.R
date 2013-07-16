QUERY_CACHE = paste(getwd(),'/','QUERY_CACHE_STANFORD','/',sep='')
DB_CFG_FILE = 'stanford_DB.cfg'
DEFAULT_DB  = 'pge_res'

QUERY_CACHE = paste(getwd(),'/','QUERY_CACHE_WHARTON','/',sep='')
DB_CFG_FILE = 'wharton_DB.cfg'
DEFAULT_DB  = 'pgefinal'

StanfordData = function(local=F){
  obj = list (
    CACHE_DIR    = paste(getwd(),'/','QUERY_CACHE_STANFORD','/',sep=''),
    DB_CFG_FILE  = 'stanford_DB.cfg',
    DEFAULT_DB   = 'pge_res',
    ZIP_DATA     = NULL,
    MULTI_PERSON = NULL
  )
  
  if(local) {
    obj$weatherDB        = function()    { db = 'pge_res'}
    obj$meterDB          = function()    { db = 'pge_res'}
    obj$resDB            = function()    { db = 'pge_res'}
    obj$meterTable       = function(zip) { table = paste('pge_res_60_',zip,sep='')}
    obj$weatherTable     = function()    { table = 'weather_60'}
    obj$accountTable     = function()    { table = 'pge_res_final'}
  } else {
    # DB specific naming and locations
    obj$weatherDB    = function()    { db = 'PGE_WEATHER'}
    obj$meterDB      = function()    { db = 'PGE_SAM'}
    obj$meterTable   = function(zip) { table = paste('pge_res_60_',zip,sep='')}
    obj$weatherTable = function()    { table = 'weather_60'}
    obj$accountTable = function()    { table = 'pge_res_final'}
  }

  # utility function that returns a list of all the zipcodes in the data set
  obj$getZips = function(useCache=F,forceRefresh=F) {
    print(obj$weatherDB())
    query    <- paste("select distinct zip5 from",obj$weatherTable(),'order by zip5')
    cacheFile = NULL
    if(useCache) { cacheFile='zipList.RData' }
    return(run.query(query,obj$weatherDB(),obj$DB_CFG_FILE,cacheDir=obj$CACHE_DIR,cacheFile=cacheFile,forceRefresh=forceRefresh)[[1]])
  }
  
  obj$getSPs = function(zip=NA,useCache=F,forceRefresh=F) {
    cacheFile = NULL
    if(is.na(zip)) { query <- paste("select distinct sp_id from",obj$meterTable()) }
    else {
      if(useCache) {
        cacheFile = paste('spids_',zip,'.RData',sep='') # only cache sp list for individual zips
      }
      query  <- paste("select distinct sp_id from",obj$meterTable(zip),"where zip5=",zip) 
    }
    return(run.query(query,obj$meterDB(),obj$DB_CFG_FILE,cacheDir=obj$CACHE_DIR,cacheFile=cacheFile,forceRefresh=forceRefresh)[[1]])
  }
  
  obj$getZipCounts = function(useCache=F) { 
    cacheFile = NULL
    if(useCache) { cacheFile=paste('zipCounts.RData',sep='') }
    query = 'SELECT zip5, COUNT(DISTINCT sp_id) as count FROM pge_res_final GROUP BY zip5'
    return(run.query(query,obj$meterDB(),obj$DB_CFG_FILE,cacheDir=obj$CACHE_DIR,cacheFile=cacheFile))
  }
  
  obj$getAllData = function(zip=NULL,useCache=F,forceRefresh=F) {
    cacheFile = NULL
    if(useCache) { cacheFile=paste('meterData_',zip,'.RData',sep='') }
    query = paste(
      'SELECT 
      sp_id, zip5, DATE,
      hkw1, hkw2, hkw3, hkw4, hkw5, hkw6, hkw7, hkw8, hkw9, hkw10,hkw11,hkw12,
      hkw13,hkw14,hkw15,hkw16,hkw17,hkw18,hkw19,hkw20,hkw21,hkw22,hkw23,hkw24 
      FROM',obj$meterTable(zip),'ORDER BY sp_id, DATE')
    zipData = run.query(query,obj$meterDB(),obj$DB_CFG_FILE,cacheDir=obj$CACHE_DIR,cacheFile=cacheFile,forceRefresh=forceRefresh)
    return(zipData)
  }
  
  obj$getSPData = function(sp_id,zip) {
    query = paste(
      'SELECT 
          sp_id,zip5,DATE,
        hkw1, hkw2, hkw3, hkw4, hkw5, hkw6, hkw7, hkw8, hkw9, hkw10,hkw11,hkw12,
        hkw13,hkw14,hkw15,hkw16,hkw17,hkw18,hkw19,hkw20,hkw21,hkw22,hkw23,hkw24 
        FROM',obj$meterTable(zip),'WHERE sp_id =',sp_id,'ORDER BY DATE')
    data = run.query(query,obj$meterDB(),obj$DB_CFG_FILE)
  }
  
  obj$getWeatherData = function(zip,useCache=F,forceRefresh=F) {
    cacheFile = NULL
    if(useCache) { cacheFile=paste('weather_',zip,'.RData',sep='') }
    query = paste(
      'SELECT `date`, TemperatureF, Pressure, DewpointF, HourlyPrecip, WindSpeed
    FROM',obj$weatherTable(),'where zip5 =',zip,'ORDER BY DATE')
    data = run.query(query,obj$weatherDB(),DATA_SOURCE$DB_CFG_FILE,cacheDir=obj$CACHE_DIR,cacheFile=cacheFile,forceRefresh=forceRefresh)
  }
  
  # return the available summary information for every zipcode
  # including sp_id count, climate zone, weather station, and income stats
  
  obj$getZipData = function(zip=NULL,useCache=F) {
    if(is.null(obj$ZIP_DATA)) {
      query = paste(
        'SELECT zip5, COUNT(DISTINCT sp_id), cecclmzn, climate, GCOUNTY, WTHRSTN,
        median_income, median_income_quantiles
        FROM', obj$accountTable(), 'GROUP BY zip5')
      # has to go into the global env to persist as a 'cache'
      cacheFile=NULL
      if(useCache) { cacheFile='zipData.RData' }
      #assign('ZIP_DATA', run.query(query,obj$resDB(),cacheFile=cacheFile), envir = .GlobalEnv)
      obj$ZIP_DATA = run.query(query,obj$resDB(),obj$DB_CFG_FILE,cacheDir=obj$CACHE_DIR,cacheFile=cacheFile)
    }
    #else { print('Using zip data cache') }
    if(is.null(zip)) { out = obj$ZIP_DATA                       }
    else             { out = obj$ZIP_DATA[obj$ZIP_DATA$zip5 == zip,]}
    return(out)
  }
  
  obj$getMultPersonSPs = function(useCache=F,forceRefresh=F) {
    if(is.null(obj$MULTI_PERSON)) {
      cacheFile = NULL
      if(useCache) { cacheFile=paste('multiPersonSPs.RData',sep='') }
      query = paste('SELECT sp_ID, COUNT(per_id) FROM pge_res_final GROUP BY sp_id HAVING COUNT(per_id) > 1')
      obj$MULTI_PERSON = run.query(query,obj$meterDB(),obj$DB_CFG_FILE,cacheDir=obj$CACHE_DIR,cacheFile=cacheFile,forceRefresh=forceRefresh)
    }
    return(obj$MULTI_PERSON)
  }
  
  class(obj) = "StanfordData"
  return(obj)
}

