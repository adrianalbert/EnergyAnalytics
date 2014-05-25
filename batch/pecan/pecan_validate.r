# pecan_validate.r
#
# Correlate ground truth end-uses with model results. 
#
# Adrian Albert
# Last modified: June 2014.
# ---------------------------------------------------------

rm(list = ls())
options(error = recover)

# __________________________________________________
# Initializations...

library('ggplot2')

DATA_PATH = '~/energy-data/pecan_street/models/'
DUMP_PATH = '~/energy-data/pecan_street/models/'
PLOT_PATH = '~/Dropbox/OccupancyStates/plots/pecan-street/validation/'

# list all data files by year/uid
files.decode = list.files(path=DATA_PATH, pattern = '*_decoded*', full.names = T, recursive = T)
files.interp = list.files(path=DUMP_PATH, pattern = '*_interpreted*', full.names = T, recursive = T)
user_info  = lapply(files.decode, function(x) {
  tmp = strsplit(x, '/')[[1]]
  ret = data.frame(res = as.character(tmp[length(tmp)-1]),                   
                   yr  = as.character(tmp[length(tmp)-3]),
                   uid = as.character(tmp[length(tmp)-2]))
  rownames(ret) = NULL
  return(ret)
})
user_info = do.call('rbind', user_info)
user_info$yr = as.character(user_info$yr)
user_info$uid= as.character(user_info$uid)
user_info$res= as.character(user_info$res)

# __________________________________________________
# Function to load data and compute metrics

file.d = files.decode[1]
file.i = files.interp[1]

compute_metrics = function(data)


