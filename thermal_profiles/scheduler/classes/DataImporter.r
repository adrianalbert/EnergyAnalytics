# #########################################################################
# DataImporter.r
# -------------------------
#
# Read in model result data. 
# 
# Adrian Albert
# Last modified: October 2013.
# #########################################################################

library('methods')
library('timeDate')
library('zoo')
library('lubridate')
require('lmtest')
library('pls')
library('reshape')

# clean-up previous definitions of methods for class Person
removeClass('DataImporter')

# ________________________
# Class definition

setClass(
  Class = "DataImporter",
  representation = representation(
    STATES_PAR     = "data.frame",          # state parameters
    TRANSITION     = "data.frame",          # transition parameters
    STATES_SES     = "data.frame"           # seasonal breakdown of states
  )
)

# function to read in data from files
read_files = function(path, pattern) {
  cur_files  = list.files(path = path, pattern = pattern, full.names = T, recursive = T)
  res        = data.frame()
  for (f in cur_files) {
    if (file.info(f)$size < 10) next
    tmp = read.csv(f)
    if (nrow(res) == 0) res = tmp else res = rbind(res, tmp)
  }
  res <- res[,colSums(is.na(res))<nrow(res)]
  res$X = NULL
  
  return(list(data = res, noChunks = length(cur_files)))
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
