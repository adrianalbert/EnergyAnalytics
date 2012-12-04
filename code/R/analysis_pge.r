# analysis_pge.r
# 
# Typical OLS analysis for PGE data:
# - regression residuals: argue non-Gaussian
# - ACF: argue non-iid
# - different models for different time-of-day/holiday/week-weekend: argue different states
# - two users: argue different responses
#
# HMM analysis (w/ covariates):
# - 
# Adrian Albert
# Last modified: November 2012.

# -----------------------------------------------
# Define constants, load libraries and functions
# -----------------------------------------------

rm(list=ls())

source('~/Dropbox/ControlPatterns/code/R/utils/timing.r')
source('~/Dropbox/ControlPatterns/code/R/utils/sql_utils.r')
source('~/Dropbox/ControlPatterns/code/R/utils/acf_ggplot.r')
source('~/Dropbox/ControlPatterns/code/R/Person.r')

N_USERS     = 10     # send batches of 50 users to each core
NOBS_THRESH = 365   # discard users with less than a year at a given premise
plots_path  = '~/Dropbox/ControlPatterns/plots/'
save_path   = '~/Dropbox/ControlPatterns/fits/'

library(utils)

# ------------------------------------------------
# Load identification information on users
# ------------------------------------------------

# SELECT * from pge_res_final3_unique WHERE total_duration >= 365 
#     INTO OUTFILE 'table.csv'
#     FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"'
#     LINES TERMINATED BY '\n';

# get a list of unique person ids
pers_info     = run.query(paste("select PER_ID, SP_ID, date from pge_res_final3_unique WHERE total_duration >", NOBS_THRESH), db = 'pge_res')
person_col    = pers_info$PER_ID
person_ids    = sort(unique(person_col))

# how many unique (PER_ID,SP_ID) tuples where people have lived for at least 1 year?
tuples        = paste(pers_info$PER_ID, pers_info$SP_ID, sep=',')
tuples_unq    = unique(tuples)
length(tuples_unq)

# divide up data into chunks
ids_vec    <- seq_along(tuples_unq)
chunks_vec <- split(tuples_unq, ceiling(ids_vec/N_USERS))

# ---------------------------------------------------------
# Wrapper to perform analysis on a given person-sp_id data
# ---------------------------------------------------------

# define some covariates 
hourly_vars   = paste('Hour.Of.Day', 0:23, sep='.')
monthly_vars  = paste('Month', 1:12, sep = '.')
weekly_vars   = paste('Day.Of.Week', c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun'), sep='.')
wthr_vars_all = c('TemperatureF', 'WindSpeed', 'Humidity', 'HourlyPrecip', 'SolarRadiation')
trend_vars    = c('Trend.8', 'Trend.24')
holiday_vars  = c('Is.Holiday.0', 'Is.Holiday.1')

source('~/Dropbox/ControlPatterns/code/R/Person.r')
personAnalysis = function(cur_data, cur_wthr, 
                           transitn.df = data.frame(), response.df = data.frame()) {
  
    # construct Person object
    user      = new(Class='Person', cur_data)  
    user      = addWeather(user, cur_wthr)
    
    cur_PER_ID = user@PER_ID
    cur_SP_ID  = user@SP_ID
    
    # _______________________________________
    # Perform analysis on current user
    
    # OLS analysis
    tic()
    user          = fitOLS(user)
    time_ols      = toc()
 
    # HMM analysis
    response_vars = c(wthr_vars_all,  weekly_vars, trend_vars) 
    transitn_vars = c(holiday_vars)
    tic()    
    user          = fitHMM(user, Kmin = 4, Kmax = 6,
                           response_vars = response_vars, transitn_vars = transitn_vars)
    time_hmm      = toc()
    
    # ______________________________
    # Save analysis results
    
    df        = as.data.frame(t(user@HMM$response$means))
    df[,setdiff(response_vars, names(df))] = NA
    df$Sigma  = user@HMM$response$stdev
    df$PER_ID = user@PER_ID
    df$SP_ID  = user@SP_ID
    df$State  = 1:user@HMM$nStates
    if (nrow(response.df)>0) response.df = rbind(response.df, df) else response.df = df
    df = as.data.frame(t(user@HMM$transition))
    df$PER_ID = user@PER_ID
    df$SP_ID  = user@SP_ID    
    if (nrow(transitn.df)>0) 
      transitn.df = rbind(transitn.df, df) else transitn.df = df

    # _______________________________
    # Produce analysis plots
    
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
    grid.newpage() 
    pushViewport(viewport(layout = grid.layout(2, 1))) 
    print(p1, vp = vplayout(1,1))
    print(p2, vp = vplayout(2,1))
    dev.off()
    
    # clear used variables
    rm(list = c('user'))
    
    return(list(transition = transitn.df, response = response.df))
}

# ------------------------------------------
# Wrapper to send computation to each core
# ------------------------------------------

analysis_wrapper <- function(chunk_id) {
  
  #sink(paste(plots_path, "log_", chunk_id,".txt", sep=''), append=TRUE)
  
  cat(paste('--------- Processing Chunk', chunk_id, '---------\n'))
  
  tic()
  per_ids   = sapply(chunks_vec[[chunk_id]], function(x) strsplit(x, ',')[[1]][1])
  sp_ids    = sapply(chunks_vec[[chunk_id]], function(x) strsplit(x, ',')[[1]][2])  
  # retrieve current zip consumption data from DB
  query     = paste("SELECT * FROM pge_res_final3_unique WHERE PER_ID IN (", paste(per_ids, collapse = ','), 
                    ') AND SP_ID IN (', paste(sp_ids, collapse = ','), 
                    ') ORDER BY date' )
  raw_data  = run.query(query, db = 'pge_res')
  zips      = unique(raw_data$ZIP5)
  
  # retrieve corresponding weather data for current zip
  query     = paste("SELECT * FROM weather_60 WHERE zip5 IN (", 
                    paste(zips, collapse = ','), ')')
  wthr_data = run.query(query, db = 'PGE_WEATHER')
  
  time_db_sec  = toc()
  
  all_users = chunks_vec[[chunk_id]]
  # perform calculations on each user in the dataset
  response.df = data.frame()
  transitn.df = data.frame()
  for (usr_info in all_users) {     

    # _________________________________
    # Select data for current person
    
    cur_PER_ID= strsplit(usr_info,',')[[1]][1]
    cur_SP_ID = strsplit(usr_info,',')[[1]][2]    
    cur_data  = subset(raw_data, PER_ID == cur_PER_ID & SP_ID == cur_SP_ID)
    zip       = unique(cur_data$ZIP5)    
    cur_wthr  = subset(wthr_data, zip5 == zip)
    
    # discard users with too few data points
    if (nrow(cur_data) < NOBS_THRESH) {
      cat(paste('--- Too few datapoints! (', nrow(cur_data), ')\n', sep=''))
      next
    }

    # perform analysis flow
    source('~/Dropbox/ControlPatterns/code/R/Person.r')
    res = try(personAnalysis(cur_data, cur_wthr, transitn.df, response.df))
    if (class(try) != 'try-error') {
      transitn.df = res$transition
      response.df = res$response
    } else {
      cat(paste('Estimation error at person', cur_PER_ID, ',', cur_SP_ID))
    }
  }
    
  # save temporary data
  write.csv(response.df, file = paste(save_path, 'response_chunk_', chunk_id, '.csv', sep=''))
  write.csv(transitn.df, file = paste(save_path, 'transitn_chunk_', chunk_id, '.csv', sep=''))
  
  # reset output connection
  #sink()
  
  rm(list = c('cur_data', 'cur_wthr'))
  return(list(response = response.df, transition = transitn.df))
}

# ------------------------------------------------
# Some preliminary analysis
# ------------------------------------------------

source('~/Dropbox/ControlPatterns/code/R/Person.r')
Rprof(filename = '~/Dropbox/ControlPatterns/fits/Rprof.out')
result = analysis_wrapper(1)
Rprof(NULL)

summaryRprof(filename = '~/Dropbox/ControlPatterns/fits/Rprof.out')

# parallel execution using multicore
library(multicore)
library(parallel)

for (chunk in 1:5){#length(chunks_vec)) {
  writeLines(c(""), paste(plots_path, "log_", chunk, ".txt", sep=''))
}
result = mclapply(1:5, FUN = analysis_wrapper, mc.silent = F, 
                   mc.preschedule = TRUE, mc.cores = 6)

