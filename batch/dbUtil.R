require(RMySQL)

# utility functions that centralize config from one database sever to the next
conf.dbCon = function(db='pge_res') {
  con = dbConnect(MySQL(), default.file=paste(baseDir,'DB_connection.cfg',sep=''), dbname=db)
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

run.query = function(query,db) {
  #print(query)
  data <- c()
  tryCatch({
    con  <- conf.dbCon(db)
    res  <- dbGetQuery(con, query)
    if(length(res)>0) data  <- res
    
  },
           error = function(e) {print(e)},
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
  return(data)
}

