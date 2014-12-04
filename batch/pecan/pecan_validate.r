# pecan_validate.r
#
# Correlate ground truth end-uses with model results. 
#
# Adrian Albert
# Last modified: June 2014.
# ---------------------------------------------------------

rm(list = ls())
options(error = recover)
options(warn=2) # warnings are treated as errors

# __________________________________________________
# Initializations...

library('ggplot2')
library('fields')
source('../../thermal_profiles/validator/metrics.r')
source('../../thermal_profiles/validator/plot_metrics.r')
source('../../utils/select_data.r')

DATA_PATH = '~/energy-data/pecan_street/usage-select/'
MODEL_PATH = '~/energy-data/pecan_street/models_new/'
PLOT_PATH = '~/Dropbox/OccupancyStates/plots/pecan-street-new/'
PRODUCE_PLOTS = T

# load user names
user_names = read.csv('~/energy-data/pecan_street/metadata/user_names_ids.csv')

# list all data files by year/uid
files.decode = list.files(path=MODEL_PATH, pattern = '*_decoded*', full.names = T, recursive = T)
files.interp = list.files(path=MODEL_PATH, pattern = '*_interpreted*', full.names = T, recursive = T)
user_info  = lapply(files.decode, function(x) {
  tmp = strsplit(x, '/')[[1]]
  ret = data.frame(res = as.character(tmp[length(tmp)-1]),                   
                   ses = as.character(tmp[length(tmp)-2]),
                   uid = as.character(tmp[length(tmp)-3]))
  rownames(ret) = NULL
  return(ret)
})
user_info            = do.call('rbind', user_info)
user_info$ses        = as.character(user_info$ses)
user_info$uid        = as.character(user_info$uid)
user_info$res        = as.character(user_info$res)

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
  tran   = lapply(unique(tran$From), function(s) {d = subset(tran, To ==s); d$From = d$To = NULL; rownames(d) = 1:nrow(d); return(d)})
  
  # check for consistency
  if (class(resp$stderr) == 'list') {
    idx.nan = which(sapply(resp$stderr, length) < 3)    
    z = c(0,0,0); names(z) = c('(Intercept)', 'TemperatureD', 'TemperatureDWinter')
    for (i in 1:length(idx.nan)) {
      tmp = resp$stderr[[idx.nan[i]]]
      z[names(tmp)] = tmp
      resp$stderr[[idx.nan[i]]] = z
    }
    resp$stderr = as.data.frame(do.call('cbind', resp$stderr))
    names(resp$stderr) = 1:ncol(resp$stderr)
  }
  
  # format model parameters for easy analysis access
  volatility = lapply(1:nStates, function(s) { z = c(mu = 0, sd = resp$stdev[s]); names(z) = c('mu', 'sd'); z})
  baseload   = lapply(1:nStates, function(s) { z = c(mu = resp$means['(Intercept)', s], sd = resp$stderr['(Intercept)', s]); names(z) = c('mu', 'sd'); z})
  response   = lapply(1:nStates, function(s) c(mu = resp$means['TemperatureD', s], sd = as.numeric(resp$stderr['TemperatureD', s])))
  bias.win   = lapply(1:nStates, function(s) {
    if ('TemperatureDWinter' %in% rownames(resp$means)) tmpVar = 'TemperatureDWinter'
    if ('TemperatureDSummer' %in% rownames(resp$means)) tmpVar = 'TemperatureDSummer'
    return(c(mu = resp$means[tmpVar, s], sd = as.numeric(resp$stderr[tmpVar, s])))
  })
  
  return(list(volatility = volatility, baseload = baseload, response = response, bias.win = bias.win,
              tran = tran, data = df))
}

# _____________________________________________________
# Function to run analysis and produce metrics & plots

run_analysis = function(data, decode_info, interp_info, PLOT_PATH = PLOT_PATH) {
  
  # assemble data
  res  = prepare_data(decode_info, data)
  df   = res$data
  df.S = select_data(df, include.weekends = FALSE, seasons = 'Summer')
  df.W = select_data(df, include.weekends = FALSE, seasons = 'Winter')
  state_names = as.character(unique(df$state))

  if (sum(is.na(df.S$AC)) >= nrow(df.S)*0.75) {
    print(" --> Not enough data for summer AC!")
    return(NULL)
  }
  if (sum(is.na(df.W$AC)) >= nrow(df.W)*0.85)  {
    print(" --> Not enough data for winter AC!")
    return(NULL)
  }
  
  if (is.null(df.S) || is.null(df.W))  {
    print(" --> Not enough data for either summer or winter")
    return(NULL)
  }

  # Compute metrics
  # -------------------------------------------------
  
  # compute state-specific ground truth linear response
  grt.S  = compute_ground_truth_response(df.S, state_names = state_names, dep.var = 'AC', indep.var = 'TemperatureD')
  grt.W  = compute_ground_truth_response(df.W, state_names = state_names, dep.var = 'AC', indep.var = 'TemperatureD')
  
  # compare state responses
  metr.S = compute_comparisons_states(grt.S, res$response)
  metr.W = compute_comparisons_states(grt.W, res$bias.win)
  
  # compute ground truth temperature-specific response curve
  temperature.grid = seq(21, 100, length.out = 80)
  a_hat.S = binned_response_curve(df.S, 
                                  bin.out = temperature.grid,
                                  dep.var = 'AC', 
                                  bin.var = 'TemperatureF', 
                                  indep.var = 'TemperatureD',
                                  nbins = 5)
  a_hat.W = binned_response_curve(df.W, 
                                  bin.out = temperature.grid,
                                  dep.var = 'AC', 
                                  bin.var = 'TemperatureF', 
                                  indep.var = 'TemperatureD',
                                  nbins = 5)
  
  # compute model and empirical state distributions
  P = compute_probability_profile(temperature.grid, res$tran)
  PeS = compute_probability_profile_empirical(df.S, 
                                              state_names = as.character(unique(df$state)),
                                              bin.out = temperature.grid,
                                              bin.var = 'TemperatureF', 
                                              nbins = 5)
  PeW = compute_probability_profile_empirical(df.W,
                                              state_names = as.character(unique(df$state)),
                                              bin.out = temperature.grid,
                                              bin.var = 'TemperatureF', 
                                              nbins = 5)

  # construct model and empirical response curves
  #a.S = model_response_curve(temperature.grid, res.S$resp, P.S)
  #a.W = model_response_curve(temperature.grid, res.W$resp, P.W)
  response.win = lapply(1:length(res$bias.win), 
                        function(i) c(res$bias.win[[i]]['mu'] + res$response[[i]]['mu'],
                                      sqrt(res$bias.win[[i]]['sd']^2 + res$response[[i]]['sd']^2)))
                  
  a.S = model_response_curve(temperature.grid, res$response, PeS$probvec)
  a.W = model_response_curve(temperature.grid, response.win, PeW$probvec)  
  
  # compute comparison metrics vs temperature (SW vs GT)
  f_vec_list = function(a) {
    lapply(1:length(a[[1]]), function(i) { z = c(a[['mu']][i], a[['sd']][i]); names(z) = c("mu", "sd"); z})
  }
  metr.temp.S  = compute_comparisons_states(f_vec_list(a_hat.S), f_vec_list(a.S))

  # reference profile to compute own error metrics
  ref.prof = lapply(1:length(a_hat.S[[1]]), function(i) c(mu = 1e-10, sd = 1e-10))
  metr.temp.S.own  = compute_comparisons_states(f_vec_list(a.S), ref.prof)
  metr.temp.Shat.own = compute_comparisons_states(f_vec_list(a_hat.S), ref.prof)
  
  # assemble return metrics: effective response by temperature
  mat= matrix(0, nrow = 8, ncol = 80); colnames(mat) = 21:100; 
  mat[1,as.character(a_hat.S$var)] = a_hat.S$mu; mat[2,as.character(a_hat.S$var)] = a_hat.S$sd
  mat[3,as.character(a_hat.W$var)] = a_hat.W$mu; mat[4,as.character(a_hat.W$var)] = a_hat.W$sd
  mat[5,as.character(a.S$var)] = a.S$mu; mat[6,as.character(a.S$var)] = a.S$sd
  mat[7,as.character(a.W$var)] = a.W$mu; mat[8,as.character(a.W$var)] = a.W$sd
  mat = as.data.frame(mat)
  mat = cbind(UID = user_ids[i], resolution = user_res[r],
              metric = rep(c('a_hat_S', 'a_hat_W', 'a_S', 'a_W'), each = 2),
              quantity = rep(c('mu', 'sd'), 4), mat)
  
  # assemble return metrics: per-state response
  r1 = cbind(metric = 'grt_S', state = 1:length(grt.S), as.data.frame(do.call('rbind', grt.S)))
  r2 = cbind(metric = 'grt_W', state = 1:length(grt.W), as.data.frame(do.call('rbind', grt.W)))
  r3 = cbind(metric = 'r_S', state = 1:length(res$response), as.data.frame(do.call('rbind', res$response)))  
  rstate = cbind(UID = user_ids[i], resolution = user_res[r], rbind(r1, r2, r3))
        
  # assemble return metrics: probabilities
  # summer magnitude
  pmat = matrix(0, nrow = nrow(t(metr.temp.S$P.magnitude)), ncol = 80); colnames(pmat) = 21:100; 
  pmat[,as.character(temperature.grid)] = t(metr.temp.S$P.magnitude)  
  p3 = as.data.frame(pmat); 
  p3 = cbind(UID = user_ids[i], resolution = user_res[r], metric = 'P_magnitude_S', eta = metr.temp.S$support.mag, p3)
  # summer relative
  pmat = matrix(0, nrow = nrow(t(metr.temp.S$P.rel.abs)), ncol = 80); colnames(pmat) = 21:100; 
  pmat[,as.character(temperature.grid)] = t(metr.temp.S$P.rel.abs)  
  p4 = as.data.frame(pmat); 
  p4 = cbind(UID = user_ids[i], resolution = user_res[r], metric = 'P_rel.abs_S', eta = metr.temp.S$support.rab, p4)
  # self error magnitude - summer
  pmat = matrix(0, nrow = nrow(t(metr.temp.S.own$P.magnitude)), ncol = 80); colnames(pmat) = 21:100; 
  pmat[,as.character(temperature.grid)] = t(metr.temp.S.own$P.magnitude)  
  p6 = as.data.frame(pmat); 
  p6 = cbind(UID = user_ids[i], resolution = user_res[r], metric = 'P_magnitude_S_own', eta = metr.temp.S.own$support.mag, p6)
  # self error magnitude - ground truth
  pmat = matrix(0, nrow = nrow(t(metr.temp.Shat.own$P.magnitude)), ncol = 80); colnames(pmat) = 21:100; 
  pmat[,as.character(temperature.grid)] = t(metr.temp.S.own$P.magnitude)  
  p7 = as.data.frame(pmat); 
  p7 = cbind(UID = user_ids[i], resolution = user_res[r], metric = 'P_magnitude_Shat_own', eta = metr.temp.Shat.own$support.mag, p7)
  
  rprob = rbind(p3, p4, p6, p7)
  metrics = list(response_temp = mat, response_state = rstate, response_prob = rprob)
      
  # Produce plots of individual performance metrics
  # -------------------------------------------------

  if (!PRODUCE_PLOTS) return(metrics)
  
  # plot density of temperature/overlap
  pdf(file=paste(PLOT_PATH, 'weather_season.pdf',sep='/'),width=9.3,height=4)
    df = data.frame(TemperatureF = c(df.W$TemperatureF, df.S$TemperatureF), Season = c(rep('Winter', nrow(df.W)), rep('Summer', nrow(df.S))))
    p = plot_overlapping_densities(df, var = 'TemperatureF', by = 'Season',
                                   title = paste('Weather for user', user_ids[i]))
    print(p)
  dev.off()

  # compare state-based response
  pdf(file=paste(PLOT_PATH, 'state_response_cmp.pdf',sep='/'),width=11,height=4)
    ll = grt.S; ll[(length(grt.S)+1):(length(grt.S)+length(res$response))] = res$response  
    ll[(length(ll)+1):(length(ll)+length(response.win))] = response.win
    p = plot_gaussian_distr(ll, 
                            state = rep(interp_info$regime.type, 3), 
                            type = rep(c('Ground', 'De-biased', 'Model'), each=length(grt.S)),
                            title = paste('Compare state responses for user', username))
    print(p)
  dev.off()

  # compare state-based response
  pdf(file=paste(PLOT_PATH, 'state_response_err_mag.pdf',sep='/'),width=9.3,height=4)
    p = plot_error_state(metr.S$P.magnitude, metr.S$support.mag, 
                         title = bquote(P(delta<eta) ~ 'for user' ~ .(username)))
    print(p)
  dev.off()
  
  # compare state-based response
  pdf(file=paste(PLOT_PATH, 'state_response_err_rel.pdf',sep='/'),width=9.3,height=4)
    rownames(metr.S$P.relative) = interp_info$regime.type
    p = plot_error_state(metr.S$P.relative, metr.S$support.rel, xlab = 'Relative Error [%/100]',
                         title = bquote(P(epsilon<eta) ~ 'for user' ~ .(username)))
    print(p)
  dev.off()
  
  pdf(file=paste(PLOT_PATH, 'state_response_err_rel_abs.pdf',sep='/'),width=9.3,height=4)
    rownames(metr.S$P.relative) = interp_info$regime.type
    p = plot_error_state(metr.S$P.rel.abs, metr.S$support.rab, xlab = 'Absolute Relative Error [%/100]',
                         title = bquote(P(abs(epsilon)<eta) ~ 'for user' ~ .(username)))
    print(p)
  dev.off()
  
  # plot state probability by temperature : empirical
  pdf(file=paste(PLOT_PATH, 'prob_profile_empiric_temp.pdf',sep='/'),width=9.3,height=4)
    colnames(PeS$probvec) = interp_info$regime.type
    p = plot_prob_profile(PeS$probvec, var = data.frame(TemperatureF = temperature.grid), 
                          title = paste('Empirical thermal regimes distribution for user', username))
    print(p)
  dev.off()  
  
  # plot state probability by temperature : model-based
  pdf(file=paste(PLOT_PATH, 'prob_profile_model_temp.pdf',sep='/'),width=9.3,height=4)
    p = plot_prob_profile(P, var = data.frame(TemperatureF = temperature.grid), 
                          title = paste('Model-based thermal regimes distribution for user', username))
    print(p)
  dev.off()  
  
  pdf(file=paste(PLOT_PATH, 'resp_curve_temp_cmp.pdf',sep='/'),width=9.3,height=4)
    p = plot_noisy_curve(list(Model = a.S, GroundTruth = a_hat.S, Winter = a.W), 
                     title = paste('Response curves for user', username))
    print(p)
  dev.off()

  pdf(file=paste(PLOT_PATH, 'resp_curve_debiased_prob_rel_abs.pdf',sep='/'),width=9.3,height=4)
    print(image.plot(temperature.grid, metr.temp.S$support.rab, metr.temp.S$P.rel.abs,
                     xlab = 'Temperature [deg F]', 
                     ylab = expression(eta), 
                     main = bquote(P(abs(epsilon)<eta) ~ '(Summer) for user' ~ .(username))))
  dev.off()
    
  pdf(file=paste(PLOT_PATH, 'resp_curve_debiased_prob_mag.pdf',sep='/'),width=9.3,height=4)
    print(image.plot(temperature.grid, metr.temp.S$support.mag, metr.temp.S$P.magnitude,
                     xlab = 'Temperature [deg F]', 
                     ylab = expression(delta), 
                     main = bquote(P(delta<eta) ~ '(Summer, de-biased) for user' ~ .(username))))
  dev.off()
  
  return(metrics)
}

# __________________________________________________
# Perform analysis across user population 

user_ids = unique(user_info$uid)
user_res = unique(user_info$res)
# user_ids = c('1801', '9922')
user_res = c('60min')
metrics = list()
for (i in 1:length(user_ids)) {       # user IDs
  for (r in 1:length(user_res)) {     # resolution levels
    
    # user name
    username = as.character(user_names$name[which(user_names$ID == user_ids[i])])
    cat(paste('Processing user', username, ' uid:', i, '/', length(user_ids), '; res:', user_res[r], '\n'))
    
    # assemble data
    f   = which(user_info$res == user_res[r] & user_info$uid == user_ids[i] & user_info$ses == 'Summer-Winter')
    
    if (length(f)==0) {
      cat(paste('user', user_ids[i], 'does not have both seasons!\n'))
      next
    }
    
    file = files.decode[f]; load(file); decode_info = data
    file.interp = files.interp[f]; load(file.interp); interp_info = data
    file.r = file.path(DATA_PATH, user_res[r], paste(user_ids[i], '.csv', sep=''))
    data   = read.csv(file.r)

    # create directory for current plots
    cur_path = file.path(PLOT_PATH, user_ids[i], user_res[r])
    if (!file.exists(cur_path)) dir.create(cur_path)
    
    # run analysis for current user
    metrics[[paste(user_ids[i], user_res[r], sep='-')]] = run_analysis(data, decode_info, interp_info, PLOT_PATH = cur_path)            
  }  
}

# save metrics to disc
save(file = file.path(MODEL_PATH, 'metrics.RData'), list = c('metrics'))
