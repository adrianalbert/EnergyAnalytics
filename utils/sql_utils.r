# sql_utils.r
# 
# Utility functions for accessing data off MySQL database.
# 
# Adrian Albert
#
# Last modified: November 2012.
# --------------------------------------------------------

# Functions to perform query on database
# -------------------------------------

require(RMySQL)

# open MySQL connection
db.conn = function(dbname = 'pge_res', verbose=F, user = 'adrian', password = 'xmenxmen') {
  if (verbose) cat('Connecting to PGE weather database...')
  db.con <- dbConnect(dbDriver("MySQL"), 
                      host   = "sssl-cyclops.stanford.edu",
                      user   = user, password = password, 	
                      dbname = dbname)
  return (db.con)
}

# Test: 
# db_cons = open.db.conn('PGE_SAM', user = 'adalbert', password = 'adrian')

# perform query
run.query = function(query,db='pge_res', verbose=FALSE, user = 'adrian', password = 'xmenxmen') {
  if (verbose) cat(paste('Performing query:',query))
  data <- c()
  tryCatch({
    con  <- open.db.conn(db, verbose = verbose, user = user, password = password)
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

# # Test:
# raw_data = run.query("select * from pge_res_final3_unique LIMIT 0,1000")
# raw_data = subset(raw_data, PER_ID == 8420562867)
