# personAnalysis.r
#
# Analysis flow for one user.
#
# Adrian Albert
# Last modified: December 2012.

# define some covariates 
hourly_vars   = paste('Hour.Of.Day', 0:23, sep='.')
monthly_vars  = paste('Month', 1:12, sep = '.')
weekly_vars   = paste('Day.Of.Week', c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun'), sep='.')
wthr_vars_all = c('TemperatureF', 'WindSpeed', 'Humidity', 'HourlyPrecip', 'SolarRadiation')
trend_vars    = c('Trend.8', 'Trend.24')
holiday_vars  = c('Is.Holiday.0', 'Is.Holiday.1')
all_covars    = c(hourly_vars, monthly_vars, weekly_vars, wthr_vars_all, trend_vars, holiday_vars)

# ---------------------------------------------------------
# Wrapper to perform analysis on a given person-sp_id data
# ---------------------------------------------------------

library('grid')
personAnalysis = function(cur_data, cur_wthr, verbose = F, ols.coefs = data.frame(), plots = T, 
                          plots_path = './plots/',
                          transitn.df = data.frame(), response.df = data.frame()) {
  
    # construct Person object
    user      = new(Class='Person', cur_data, log = T)  
    user      = addWeather(user, cur_wthr)
    
    cur_PER_ID = user@PER_ID
    cur_SP_ID  = user@SP_ID
    
    # _______________________________________
    # Perform analysis on current user
    
    # OLS analysis
    user          = fitOLS(user, verbose = verbose)
 
    # HMM analysis
    response_vars = c(wthr_vars_all,  hourly_vars, trend_vars, weekly_vars) 
    transitn_vars = c()
    user          = fitHMM(user, Kmin = 5, Kmax = 5, constrMC = NULL, verbose = verbose, 
                           response_vars = response_vars, transitn_vars = transitn_vars)

    # ______________________________
    # Save analysis results: OLS
    
    df = as.data.frame(t(user@OLS$coef.signif[,1]))
    df[,setdiff(all_covars, names(df))] = NA
    df = df[,all_covars]
    df$R2 = user@OLS$fit.summary$adj.r.squared
    df$PER_ID = user@PER_ID
    df$SP_ID  = user@SP_ID
    df$ZIPCODE= user@ZIPCODE
    if (nrow(ols.coefs)>0) ols.coefs = rbind(ols.coefs, df) else ols.coefs = df
    
    # ______________________________
    # Save analysis results: HMM
    
    df        = as.data.frame(t(user@HMM$response$means))
    df[,setdiff(response_vars, names(df))] = NA
    df$Sigma  = user@HMM$response$stdev
    df$PER_ID = user@PER_ID
    df$SP_ID  = user@SP_ID
    df$ZIPCODE= user@ZIPCODE
    df$State  = 1:user@HMM$nStates
    if (nrow(response.df)>0) response.df = rbind(response.df, df) else response.df = df
    df = as.data.frame(t(user@HMM$transition))
    df$Transit= rownames(df)
    rownames(df) = NULL
    df$PER_ID = user@PER_ID
    df$SP_ID  = user@SP_ID    
    df$ZIPCODE= user@ZIPCODE
    if (nrow(transitn.df)>0) 
      transitn.df = rbind(transitn.df, df) else transitn.df = df

    # _______________________________
    # Produce analysis plots
    
    if (plots) {
      interval = c(min(user@timestamps),  
                   min(user@timestamps) + 3600 * 24 * 7) + 3600 * 24 * 5    
  
      # OLS fit
      png(paste(plots_path, paste(cur_PER_ID, cur_SP_ID, sep='_'), '_OLS_fit.png', sep=''), width=1200, height=600)
      plot(user, type='OLS-fit', interval=interval)
      dev.off()
      
      # OLS residuals
      png(paste(plots_path, paste(cur_PER_ID, cur_SP_ID, sep='_'), '_OLS_resid.png', sep=''), width=1200, height=600)
      plot(user, type='OLS-res', interval=interval)
      dev.off()
  
      # HMM fit
      p1 = plot(user, type='HMM-ts', interval = interval)
      png(paste(plots_path, paste(cur_PER_ID, cur_SP_ID, sep='_'), '_HMM_fit.png', sep=''), width=1000, height=400)
      print(p1)
      dev.off()
      
      # HMM analysis
      p2 = plot(user, type='HMM-MC')
      p3 = plot(user, type='HMM-ci')    
      p4 = plot(user, type='HMM-acf')    
      png(paste(plots_path, paste(cur_PER_ID, cur_SP_ID, sep='_'), '_HMM_analysis.png', sep=''), width=1000, height=600)
      grid.newpage() 
      pushViewport(viewport(layout = grid.layout(2, 2))) 
      vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y) 
      print(p4, vp = vplayout(1, 1:2))
      print(p3, vp = vplayout(2,1))
      print(p2, vp = vplayout(2,2))
      dev.off()
  
      # HMM residuals
      p1 = plot(user, type = 'HMM-res')
      p2 = plot(user, type = 'HMM-qq')
      png(paste(plots_path, paste(cur_PER_ID, cur_SP_ID, sep='_'), '_HMM_res.png', sep=''), width=1200, height=500)
      print(p1)
      grid.newpage() 
      pushViewport(viewport(layout = grid.layout(2, 1))) 
      print(p1, vp = vplayout(1,1))
      print(p2, vp = vplayout(2,1))
      dev.off()
    }
    
    # clear used variables
    rm(list = c('user'))
    
    return(list(transition = transitn.df, response = response.df, ols.coefs = ols.coefs))
}

