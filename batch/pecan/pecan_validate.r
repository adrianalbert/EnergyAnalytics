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
source('../../thermal_profiles/validator/')


DATA_PATH = '~/energy-data/pecan_street/usage-select/'
DUMP_PATH = '~/energy-data/pecan_street/models/'
PLOT_PATH = '~/Dropbox/OccupancyStates/plots/pecan-street/validation/'

# list all data files by year/uid
files.decode = list.files(path=DUMP_PATH, pattern = '*_decoded*', full.names = T, recursive = T)
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

f = 1
file.d = files.decode[f]; load(file.d); decode_info = data
file.i = files.interp[f]; load(file.i); interp_info = data
file.r = file.path(DATA_PATH, user_info[f, 'res'], user_info[f, 'yr'], paste(user_info[f, 'uid'], '.csv', sep=''))
data   = read.csv(file.r)

# prepare data for analysis
prepare_data = function(decode_info, interp_info, data) {
  
  # access data
  nTrain = length(decode_info$state)
  nStates= decode_info$nStates
  df     = data[1:nTrain, ]; df$state = decode_info$states; 
  df$TemperatureD = df$TemperatureF - 65
  resp   = decode_info$response   
  tran   = decode_info$transition
  tran   = lapply(unique(tran$From), function(s) {d = subset(tran, From ==s); d$From = d$To = NULL; rownames(d) = 1:nrow(d); return(d)})
  
  return(list(resp = resp, tran = tran, data = df))
}
  
# function to run analysis at different levels of data resolution


res = prepare_data(decode_info, interp_info, data)
df = res$data
