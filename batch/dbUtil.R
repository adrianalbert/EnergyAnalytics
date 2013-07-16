require(RMySQL)

QUERY_CACHE = paste(getwd(),'/','QUERY_CACHE','/',sep='')
DB_CFG_FILE = 'stanford_DB.cfg'
DEFAULT_DB  = 'pge_res'

# utility functions that centralize config from one database sever to the next
conf.dbCon = function(db=DEFAULT_DB,cfgFile=DB_CFG_FILE) {
  con = dbConnect(MySQL(), default.file=paste(getwd(),'/',cfgFile,sep=''), dbname=db)
  return(con)
}

# utility fn to clear all active variables - 
#leaves .varName vars behind to eliminate these, use ls(all.names=TRUE)
clear = function() { rm(list=ls(),envir=baseenv()) }

# utility fn to close/clear all active db connections
clearCons = function() {
  all_cons <- dbListConnections(MySQL())
  for(con in all_cons) {
    res <- dbListResults(con)
    if(length(res)>0) { 
      dbClearResult(res[[1]]) 
      rm(res)
    }
    dbDisconnect(con) 
  }
}

# use this when too many cons are open
showCons = function() {
  all_cons <- dbListConnections(MySQL())
  print(dim(all_cons))
  print(all_cons)
  s = dbGetQuery(all_cons[[1]], "show processlist") 
  print(s)
}

run.query = function(query,db,cfgFile,cacheDir=NULL,cacheFile=NULL,forceRefresh=F) {
  QUERY_RESULT <- c()
  cachePath = NULL
  if(! is.null(cacheFile)) {
    # try to load from disk
    dir.create(file.path(cacheDir),showWarnings=FALSE)
    cachePath = paste(cacheDir,cacheFile,sep='')
    if(file.exists(cachePath)) {
      print(paste('Data cache found. Loading data from',cacheFile))
      load(cachePath) # this should load into the variable QUERY_RESULT
    }
  }
  if(length(QUERY_RESULT) == 0 | forceRefresh) { # skip the DB stuff if it has been loaded from disk
    print(query)
    con = NULL
    tryCatch({
      con  <- conf.dbCon(db,cfgFile)
      res  <- dbGetQuery(con, query)
      if(length(res)>0) QUERY_RESULT  <- res
    },
    error = function(e) {print(paste('Error in run.query:',e))},
    finally = {
      # close the results set if necessary
      resultSet <- dbListResults(con)
      if(length(resultSet)>0) { 
        dbClearResult(resultSet[[1]])
        rm(resultSet)
      }
      dbDisconnect(con)
      rm(con)
    } )
    if(! is.null(cachePath)) {          # save results to disk cache for future use
      save(QUERY_RESULT,file=cachePath)
    }
  }
  return(QUERY_RESULT)
}

