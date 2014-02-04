# #########################################################################
# DataImporter.r
# -------------------------
#
# Read in model result data. 
# 
# Adrian Albert
# Last modified: February 2014.
# #########################################################################

library('methods')
library('timeDate')
library('zoo')
library('lubridate')
require('lmtest')
library('pls')
library('reshape')

# clean-up previous definitions of methods for class 
removeClass('DataImporter')

# ________________________
# Class definition

setClass(
  Class = "DataImporter",
  representation = representation(
    STATES_PAR     = "data.frame",          # state parameters
    TRANSITION     = "data.frame"           # transition parameters
  )
)

# function to read in data from RData files
read_rdata_files = function(path) {
  
  # read decoder data
  dec_files  = list.files(path = path, pattern = '*decoded*', full.names = T, recursive = T)
  int_files  = list.files(path = path, pattern = '*interpreted*', full.names = T, recursive = T)
  cur_uids   = list.files(path = path, pattern = '*interpreted*', full.names = F, recursive = T)  
  cur_uids   = sapply(strsplit(cur_uids, '_'), function(l) as.numeric(l[[1]]))
  resp       = list()
  sder       = list()
  init       = list()
  for (i in 1:length(cur_uids)) {
    load(dec_files[i])
    resp[[1+length(resp)]] = cbind(data$response$means, UID = cur_uids[i])
    sder[[1+length(sder)]] = cbind(rbind(stdev = data$response$stdev, stderr = data$response$stderr[2,]), UID = cur_uids[i])
    load(int_files[i])
    init[[1+length(init)]] = cbind(data$benchmarks$pi0, UID = cur_uids[i])
  }
  
  # 
  
  return(list(data = res))
}
# _____________________________________________________
# Constructor method for class DataImporter

setMethod(f = "initialize", 
          signature = "DataImporter",
          definition = function(.Object, 
                                path = './', 
                                states_par = 'response',
                                transition = 'transitn',
                                states_ses = 'hmmseasn',
                                verbose = T) {
            
            if (verbose) cat('*** Reading in data ***\n')
            
            # ____________________________
            # Format state attribute data
            
            res = read_files(path, states_par)
            .Object@STATES_PAR = res[['data']]
            if (verbose) cat(paste('State parameters:\tRead in', res[['noChunks']], 'chunks\n'))
            
            # ____________________________
            # Format transition data
            
            res = read_files(path, transition)
            .Object@TRANSITION = res[['data']]
            if (verbose) cat(paste('Transition parameters:\tRead in', res[['noChunks']], 'chunks\n'))              
            # ____________________________
            # Format benchmark data
            
            res = read_files(path, states_ses)
            .Object@STATES_SES = res[['data']]
            if (verbose) cat(paste('Seasonal breakdown:\tRead in', res[['noChunks']], 'chunks\n'))              

            return(.Object)
          })


# ____________________________________________________________
# Print method for class DataFormatter: display useful stats.

setMethod('show',
          signature  = 'DataImporter',
          definition = function(object){
            cat('*** DataImporter Object ***\n')
            
            # basic info
            cat(sprintf('No. Users:   %d\n', length(unique(object@STATES_PAR$UID))))            
            
            cat('*** END: DataImporter Object ***\n')
          })
