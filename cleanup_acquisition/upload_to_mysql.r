# --------------------------------------
# Script to upload data in csv files to 
# MySQL tables.
# 
# Adrian Albert
# October 2012.
# --------------------------------------

rm(list=ls())

# Read files and current tables in database
# ------------------------------------------

# Read file names from data folder
path      = "~/ControlPatterns/data/PGE_WEATHER_RESIDENTIAL/"
all.files = list.files(path = path, pattern = '[0-9]{4}', all.files = FALSE,
                       full.names = FALSE, recursive = FALSE,
                       ignore.case = FALSE, include.dirs = FALSE)
tabs.new  = sapply(all.files, function(file) {
  zip = strsplit(file, '\\.')[[1]][1]
  return(paste('ZIP_',zip,sep=''))
})

# open a connection to MySQL database
library(RMySQL) # will load DBI as well
con <- dbConnect(dbDriver("MySQL"), 
                 host = "sssl-cyclops.stanford.edu",
                 user = "adrian", password = "mypassword", 
                 dbname = "PGE_WEATHER")

# list the tables in the database
tabs.old = dbListTables(con)

# only upload new tables to database
upload.tabs = setdiff(tabs.new, tabs.old)

# perform queries
for (tab in upload.tabs) {
  zip  = strsplit(tab, '_')[[1]][1] 
  file = paste(path, zip, '.csv', sep='')
  cur.table = read.csv(file)
  dbWriteTable(con, 
               tab, 
               cur.table, 
               overwrite = TRUE)  
}
dbDisconnect(con)                     

