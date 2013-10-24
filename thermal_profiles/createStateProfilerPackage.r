# createStateProfilerPackage.r
#
# Create R package "stateProfiler"
# 
# Adrian Albert
# Last modified: October 2013.
# -----------------------------------------------------------------------
# 
rm(list = ls())
files = c(list.files(path = "./classes", full.names=T), 
          list.files(path = "./visualization", full.names=T), 
          './depmixfit-modified.R', 'viterbi_states.R')
package.skeleton(name = 'stateProfiler', code_files = files)