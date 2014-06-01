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
library('fields')
source('../../thermal_profiles/validator/metrics.r')
source('../../thermal_profiles/validator/plot_metrics.r')

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
# Function to load data and prepare analysis

# prepare data for analysis
prepare_data = function(decode_info, data) {
  
  # access data
  nTrain = length(decode_info$state)
  nStates= decode_info$nStates
  df     = data[1:nTrain, ]; df$state = decode_info$states; 
  df$TemperatureD = df$TemperatureF - 65
  resp   = decode_info$response   
  tran   = decode_info$transition
  tran   = lapply(unique(tran$From), function(s) {d = subset(tran, From ==s); d$From = d$To = NULL; rownames(d) = 1:nrow(d); return(d)})
  
  # format model parameters for easy analysis access
  volatility = resp$stdev
  baseload   = lapply(1:nStates, function(s) c(mu = resp$means['(Intercept)', s], sd = resp$stderr['(Intercept)', s]))
  response   = lapply(1:nStates, function(s) c(mu = resp$means['TemperatureD', s], sd = as.numeric(resp$stderr['TemperatureD', s])))
  
  return(list(volatility = volatility, baseload = baseload, response = response, 
              tran = tran, data = df))
}


f = 1
file.d = files.decode[f]; load(file.d); decode_info = data
file.i = files.interp[f]; load(file.i); interp_info = data
file.r = file.path(DATA_PATH, user_info[f, 'res'], user_info[f, 'yr'], paste(user_info[f, 'uid'], '.csv', sep=''))
data   = read.csv(file.r)

res = prepare_data(decode_info, data)
df  = res$data

# __________________________________________________
# Function to compute performance metrics

source('../../thermal_profiles/validator/metrics.r')

# compute state-specific ground truth linear response
grt = compute_ground_truth_response(df, dep.var = 'nonHVAC', indep.var = 'TemperatureD')

# compare state responses
metr= compute_comparisons_states(grt, res$response)

# compute ground truth temperature-specific response curve
temperature.grid = seq(40, 100, length.out = 50)
a_hat = binned_response_curve(df, 
                             bin.out = temperature.grid,
                             dep.var = 'nonHVAC', 
                             bin.var = 'TemperatureF', 
                             indep.var = 'TemperatureD',
                             nbins = 6)

# compute model temperature curve
P = compute_probability_profile(temperature.grid, res$tran)
a = model_response_curve(temperature.grid, res$resp, res$tran)

f_vec_list = function(a) {
  lapply(1:length(a[[1]]), function(i) c(mu = a[['mu']][i], sd = a[['sd']][i]))
}
metr.temp = compute_comparisons_states(f_vec_list(a_hat), f_vec_list(a))

# __________________________________________________
# Produce plots of individual performance metrics

source('../../thermal_profiles/validator/plot_metrics.r')
ll = grt; ll[(length(grt)+1):(length(grt)+length(res$response))] = res$response
plot_gaussian_distr(ll, 
                    state = as.factor(rep(1:length(grt), 2)), 
                    type = rep(c('Ground', 'Model'), each=length(grt)))

plot_error_state(metr$P.magnitude, metr$support.mag)

plot_error_state(metr$P.relative, metr$support.rel, xlab = 'Relative Error [%/100]')

plot_error_state(metr$P.rel.abs, metr$support.rab, xlab = 'Absolute Relative Error [%/100]')

plot_noisy_curve(list(Model = a, GroundTruth = a_hat), 
                 title = 'Model and ground truth response curves')

image.plot(temperature.grid, metr.temp$support.rab, metr.temp$P.rel.abs,
           xlab = 'Temperature [deg F]', 
           ylab = expression(eta), 
           main = expression(P(abs(epsilon)<eta)))

# ____________________________________________________________
# Function to perform analysis at different resolution levels


