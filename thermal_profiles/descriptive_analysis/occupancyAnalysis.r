# occupancyAnalysis.r
#
# Analysis flow for one user.
#
# Adrian Albert
# Last modified: May 2013.

# ---------------------------------------------------------
# Wrapper to perform analysis on a given time series data.
# ---------------------------------------------------------

options(error = recover)

source('code/OccupancyStates.r')
library('grid')

occupancyAnalysis = function(cur_data, cur_wthr, UID, zip, 
                             NOBS_THRESH = 90, Kmin = 3, Kmax = 6,
                             maxit = 100, nRestarts = 2,
                             thresh.R2 = 0.8, thresh.MAPE = 0.15,
                             verbose = F, plots_path = NULL, dump_path = NULL,
                             resp.vars = c('(Intercept)', 'TemperatureF'), 
                             tran.vars = c('(Intercept)'),
                             addl.vars = c('Weekend')) {
  
  if (!is.null(dump_path)) if (!file.exists(dump_path))  dir.create(dump_path, recursive = TRUE)
  
  # construct OccupancyStates object
  user      = new(Class='OccupancyStates', cur_data, UID, zip, verbose = verbose)  
  if (ncol(cur_wthr) == 2) {    
    wthr_data = data.frame(cur_wthr[,-which(names(cur_wthr) == 'date')])
    names(wthr_data) = setdiff(names(cur_wthr), 'date')
  } else wthr_data = cur_wthr[,-which(names(cur_wthr) == 'date')]
  
  user      = addWeather(user, cur_wthr$date, wthr_data, verbose = verbose)    
  
  if (length(user@kWh) < NOBS_THRESH * 24) {
    stop(paste('---> Too few datapoints at user ', user@UID, '! (', round(length(user@kWh)/24), ')\n', sep=''))
  }
  
  # prepare data
  user      = prepareData(user, train.frac = 0.9, 
                          resp.vars = resp.vars,
                          tran.vars = tran.vars,
                          addl.vars = addl.vars)    
  if (nrow(user@data.train) == 0) stop(paste('---> No useful data for user', user@UID))
  
  # OLS analysis
  user          = fitOLS(user, verbose = verbose)
  user          = computeBreakpoint(user, break.var = c(TemperatureF=55, TemperatureF=60))    
  
  # HMM analysis
  user          = learnOccupancyStates(user, Kmin = Kmin, Kmax = Kmax, 
                                       maxit = maxit, nRestarts = nRestarts,
                                       thresh.R2 = thresh.R2, thresh.MAPE = thresh.MAPE)
  show(user)

  user          = computePredictionAccuracy(user, test.periods = 12, verbose = verbose)
  user          = computeStatsHMM(user)

  # compute useful metrics from model
  user          = computeMetricsHMM(user)
  
#  user          = pruneStates(user, thresholds = c(none = 0.10, low = 0.5, high = 2))
  user          = computeContributionsHMM(user, verbose = verbose)

  show(user)
  
  # ______________________________
  # Save analysis results
  
  res1      = dumpComputationToFile(user, path = NULL)
  res2      = dumpComputation(user, path = NULL)    
  
  # _______________________________
  # Produce analysis plots
  
  if (!is.null(plots_path)) {
    if (!file.exists(plots_path)) dir.create(plots_path, recursive = TRUE)    
    timestamps = as.POSIXct(user@timestamps)
    interval = c(min(timestamps),  
                 min(timestamps) + 3600 * 24 * 7) + 3600 * 24 * 5    

    # OLS fit
#     png(paste(plots_path, user@UID, '_OLS_fit.png', sep=''), width=1200, height=600, res = 100)
#     plot(user, type='OLS-fit', interval=interval)
#     dev.off()
    
    # OLS residuals
    png(paste(plots_path, user@UID, '_OLS_resid.png', sep=''), width=1200, height=600, res = 100)
    plot(user, type='OLS-res', interval=interval)
    dev.off()

    # HMM fit
    p1 = plot(user, type='HMM-ts', interval = interval)
    png(paste(plots_path, user@UID, '_HMM_fit.png', sep=''), width=1400, height=600, res = 100)
    print(p1)
    dev.off()

    # heatmap plots
    png(paste(plots_path, user@UID, '_HMM-consumption-heatmap.png', sep=''), width = 1400, height = 800, res = 100)
    plot(user, type = 'HMM-state-heatmap', covar = 'kWh')
    dev.off()
    
    png(paste(plots_path, user@UID, '_HMM-temperature-heatmap.png', sep=''), width = 1400, height = 800, res = 100)
    plot(user, type = 'HMM-state-heatmap', covar = 'TemperatureF')
    dev.off()
    
    png(paste(plots_path, user@UID, '_HMM-state-heatmap.png', sep=''), width = 1400, height = 800, res = 100)
    plot(user, type = 'HMM-state-heatmap', covar = 'States')
    dev.off()
    
    png(paste(plots_path, user@UID,'_HMM-state-breakdown.png', sep=''), width = 1400, height = 800, res = 100)
    print(plot(user, type = 'HMM-state-breakdown'))
    dev.off()
    
    png(paste(plots_path, user@UID, '_HMM-state-heatmap2.png', sep=''), width = 2000, height = 1000, res = 200)
    print(plot(user, type = 'HMM-state-heatmap2', covar = 'States'))
    dev.off()      
    
    # coefficient time series
    png(paste(plots_path, user@UID, '_HMM-coefs-ts.png', sep=''), width = 1400, height = 800, res = 100)
    print(plot(user, type = 'HMM-coefs-ts', interval = interval, covar = c('TemperatureF')))
    dev.off()      
      
    # HMM structure
    png(paste(plots_path, user@UID, '_HMM_MC_cov.png', sep=''), width=1200, height=800, res = 100)
    print(plot(user, type='HMM-MC-cov'))
    dev.off()
    
    # HMM PACF
    png(paste(plots_path, user@UID, '_HMM_pacf.png', sep=''), width=1000, height=600, res = 100)
    print(plot(user, type='HMM-pacf', PACF=T))
    dev.off()
      
    # HMM residuals
    png(paste(plots_path, user@UID, '_HMM_res.png', sep=''), width=1200, height=500, res = 100)
    plot(user, type = 'HMM-res')
    dev.off()
    
    # HMM covariate contributions 
    png(paste(plots_path, user@UID, '_HMM_contrib_ts.png', sep=''), width = 1400, height = 800, res = 100)
    plot(user, type = 'HMM-contrib-ts', interval=interval, covar = c('TemperatureF'))
    dev.off()

    # HMM covariate contributions (total)
    png(paste(plots_path, user@UID, '_HMM_contrib_tot.png', sep=''), width = 800, height = 400, res = 100)
    print(plot(user, type = 'HMM-contrib-tot', nrow=2))
    dev.off()      
    
    # temperature dependence
    png(paste(plots_path, user@UID, '_HMM-dep-covar.png', sep = ''), width = 1000, height = 600, res = 150)
    print(plot(user, type = 'HMM-dep-covar', covar = 'TemperatureF'))
    dev.off()      
    
    png(paste(plots_path, user@UID, 'HMM-dep-covar-sep.png', sep=''), width = 1000, height = 600, res = 150)
    print(plot(user, type = 'HMM-dep-covar', covar = 'TemperatureF', separate = T))
    dev.off()
    
    png(paste(plots_path, user@UID, '_HMM-dep-covar-orig.png', sep = ''), width = 1000, height = 600, res = 150)
    print(plot(user, type = 'HMM-dep-covar', covar = 'TemperatureF', highlight = F))
    dev.off()
    
    png(paste(plots_path, user@UID,'HMM-trans-prob.png',sep=''), width = 1600, height = 800, res = 100)
    print(plot(user, type = 'HMM-trans-prob', covar = 'TemperatureF'))
    dev.off()
    
    png(paste(plots_path, user@UID,'HMM-stationary-prob.png',sep=''), width = 1000, height = 400, res = 100)
    print(plot(user, type = 'HMM-stationary-prob', covar = 'TemperatureF'))
    dev.off()
    
    png(paste(plots_path, user@UID,'HMM-err-horiz.png',sep=''), width = 1000, height = 400, res = 100)
    print(plot(user, type = 'HMM-err-horiz'))
    dev.off()
    
    png(paste(plots_path, user@UID,'HMM-time-temp.png',sep=''), width = 1000, height = 400, res = 100)
    print(plot(user, type = 'HMM-time-temp', covar = 'TemperatureF'))
    dev.off()
    
    png(paste(plots_path, user@UID,'HMM-aggregate-season.png',sep=''), width = 1600, height = 600, res = 100)
    print(plot(user, type = 'HMM-aggregate-season'))
    dev.off()
    
  }
  
  return(list(user,res2,res1))
}

