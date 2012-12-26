# personAnalysis.r
#
# Analysis flow for one user.
#
# Adrian Albert
# Last modified: December 2012.

# ---------------------------------------------------------
# Wrapper to perform analysis on a given person-sp_id data
# ---------------------------------------------------------

# define some covariates 
hourly_vars   = paste('HourOfDay', 0:23, sep='.')
monthly_vars  = paste('Month', 1:12, sep = '.')
weekly_vars   = paste('DayOfWeek', c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun'), sep='.')
wthr_vars_all = c('TemperatureF', 'WindSpeed', 'Humidity', 'HourlyPrecip', 'SolarRadiation')
trend_vars    = c('Trend.8', 'Trend.24')
holiday_vars  = c('IsHoliday.0', 'IsHoliday.1')
all_covars    = c(hourly_vars, monthly_vars, weekly_vars, wthr_vars_all, trend_vars, holiday_vars)

library('grid')
personAnalysis = function(cur_data, cur_wthr, UID, zip, 
                          verbose = F, ols.coefs = data.frame(), plots = T,
                          plots_path = './plots/', fits_path = './fits/',
                          transitn.df = data.frame(), response.df = data.frame(), 
			  ols.comps = data.frame(), hmm.comps = data.frame()) {
      
    # construct Person object
    user      = new(Class='Person', cur_data, UID, zip, log = T, verbose = verbose)  

print(range(user@timestamps))	

    user      = addWeather(user, cur_wthr, verbose = verbose)    
    
    # _______________________________________
    # Perform analysis on current user
    
    # prepare data
    user          = prepareData(user, trends = c(8,24,24*6*30))
    if (nrow(user@data.tmp) == 0) stop(paste('---> No useful data for user', user@UID))
    
    # OLS analysis
    user          = fitOLS(user, verbose = verbose, stats = TRUE)
 
    # HMM analysis
    response_vars = c(wthr_vars_all,  hourly_vars, trend_vars, weekly_vars) 
    transitn_vars = c()
    user          = fitHMM(user, Kmin = 5, Kmax = 5, constrMC = NULL, verbose = verbose, ols_vars = T,
                           response_vars = response_vars, transitn_vars = transitn_vars)

    # compute covariate contributions
    user = computeContributionsOLS(user, verbose = verbose)
    user = computeContributionsHMM(user, verbose = verbose)

    show(user)
    
    # ______________________________
    # Save analysis results
    
    # format data for saving: coefficients
    res  = dumpComputationToFile(user, path = fits_path, dump = FALSE)
    
    if (nrow(ols.coefs)>0) ols.coefs = rbind(ols.coefs, res$OLS) else ols.coefs = res$OLS    
    if (nrow(response.df)>0) response.df = rbind(response.df, res$HMM.response) else response.df = res$HMM.response
    if (nrow(transitn.df)>0) transitn.df = rbind(transitn.df, res$HMM.transition) else transitn.df = res$HMM.transition
    if (nrow(ols.comps)>0) ols.comps = rbind(ols.comps, res$OLS.comp) else ols.comps = res$OLS.comp
    if (nrow(hmm.comps)>0) hmm.comps = rbind(hmm.comps, res$HMM.comp) else hmm.comps = res$HMM.comp

    # _______________________________
    # Produce analysis plots
    
    if (plots) {
      timestamps = as.POSIXct(user@timestamps)
      interval = c(min(timestamps),  
                   min(timestamps) + 3600 * 24 * 7) + 3600 * 24 * 5    
  
      # OLS fit
      png(paste(plots_path, user@UID, '_OLS_fit.png', sep=''), width=1200, height=600)
      plot(user, type='OLS-fit', interval=interval)
      dev.off()
      
      # OLS residuals
      png(paste(plots_path, user@UID, '_OLS_resid.png', sep=''), width=1200, height=600)
      plot(user, type='OLS-res', interval=interval)
      dev.off()
  
      # OLS covariate contributions
      png(paste(plots_path, user@UID, '_OLS_contrib_ts.png', sep=''), width = 1200, height = 500)
      print(plot(user, type = 'OLS-contrib-ts', interval=interval))
      dev.off()

      # OLS covariate contributions (total)
      png(paste(plots_path, user@UID, '_OLS_contrib_tot.png', sep=''), width = 800, height = 600)
      print(plot(user, type = 'OLS-contrib-tot'))
      dev.off()

      # HMM fit
      p1 = plot(user, type='HMM-ts', interval = interval)
      png(paste(plots_path, user@UID, '_HMM_fit.png', sep=''), width=1000, height=400)
      print(p1)
      dev.off()
      
      # HMM analysis
      p2 = plot(user, type='HMM-MC')
      p3 = plot(user, type='HMM-ci')    
      p4 = plot(user, type='HMM-acf')    
      png(paste(plots_path, user@UID, '_HMM_analysis.png', sep=''), width=1000, height=600)
      grid.newpage() 
      pushViewport(viewport(layout = grid.layout(2, 2))) 
      vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y) 
      print(p4, vp = vplayout(1, 1:2))
      print(p3, vp = vplayout(2,1))
      print(p2, vp = vplayout(2,2))
      dev.off()
      
      # MC structure
      png(paste(plots_path, user@UID, '_HMM_MC.png', sep=''), width=1000, height=600)
      plot(user, type='HMM-MC2')
      dev.off()
      
      # HMM residuals
      png(paste(plots_path, user@UID, '_HMM_res.png', sep=''), width=1200, height=500)
      plot(user, type = 'HMM-res')
      dev.off()
      
      # HMM covariate contributions 
      png(paste(plots_path, user@UID, '_HMM_contrib_ts.png', sep=''), width = 1200, height = 500)
      print(plot(user, type = 'HMM-contrib-ts', interval=interval))
      dev.off()

      # HMM covariate contributions (total)
      png(paste(plots_path, user@UID, '_HMM_contrib_tot.png', sep=''), width = 1200, height = 500)
      print(plot(user, type = 'HMM-contrib-tot', nrow=2))
      dev.off()
    }
    
    # clear used variables
    rm(list = c('user'))
    
    return(list(transition = transitn.df, response = response.df, ols.coefs = ols.coefs,
		 ols.comps = ols.comps, hmm.comps = hmm.comps))
}

