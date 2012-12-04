# fit_reg_hmm_pge.r
# 
# Fits HMM to PGE data.
#
# Adrian Albert
# Last modified: November 2012.

rm(list=ls())

# Access data in database
# -----------------------

library(RMySQL) # will load DBI as well

cat('Connecting to PGE consumption database...')
cons_db <- dbConnect(dbDriver("MySQL"), 
                 host = "sssl-cyclops.stanford.edu",
                 user = "adrian", password = "gambit", 
                 dbname = "PGE_SAM")
cat('Loading consumption tables list...')
cons_tabs = dbListTables(cons_db)

cat('Connecting to PGE weather database...')
wthr_db <- dbConnect(dbDriver("MySQL"), 
                 host = "sssl-cyclops.stanford.edu",
                 user = "adrian", password = "gambit", 
                 dbname = "PGE_WEATHER")
cat('Loading consumption tables list...')
wthr_tabs = dbListTables(wthr_db)

query  = paste("select * from ", wthr_tabs[16], "LIMIT 0, 5")
weather= dbGetQuery(wthr_db, query)

query  = paste("select * from ", cons_tabs[3])
test   = dbGetQuery(cons_db, query)

# Fit linear regression model
# ---------------------------

# Analyze OLS residuals
# ---------------------


cat('Disconnecting from PGE database...')
dbDisconnect(con)                     
