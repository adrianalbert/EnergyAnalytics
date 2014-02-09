# #########################################################################
# DataImporter.r
# -------------------------
#
# Read in model result data. 
# 
# Adrian Albert
# Last modified: February 2014.
# #########################################################################

library('methods')
library('timeDate')
library('zoo')
library('lubridate')
require('lmtest')
library('pls')
library('reshape')
library('fields')
source('../../utils/viz/plot_utils.r')

# clean-up previous definitions of methods for class 
removeClass('DataImporter')

# ________________________
# Class definition

setClass(
  Class = "DataImporter",
  representation = representation(
    RESPONSE        = "data.frame",          # state parameters
    RESP.ERR        = "data.frame",          # state parameters errors
    VOLATILITY      = "data.frame",          # state-specific volatilities
    INIT.DIS        = "data.frame",          # initial distribution over states
    TRANSITION      = "data.frame",          # transition parameters
    PERFORMANCE     = "data.frame"           # performance stats
  )
)

# function to read in data from RData files
read_rdata_files = function(path) {
  
  # read decoder data
  dec_files  = list.files(path = path, pattern = '*decoded*', full.names = T, recursive = T)
  int_files  = list.files(path = path, pattern = '*interpreted*', full.names = T, recursive = T)
  cur_uids   = list.files(path = path, pattern = '*interpreted*', full.names = F, recursive = T)  
  cur_uids   = sapply(strsplit(cur_uids, '_'), function(l) as.character(l[[1]]))
  resp       = list()
  serr       = list()
  sdev       = list()
  init       = list()
  tran       = list()
  perf       = list()  
  for (i in 1:length(cur_uids)) {
    
    load(dec_files[i])
    K = ncol(data$response$means)
    resp[[1+length(resp)]] = cbind(data$response$means, UID = cur_uids[i], variable = rownames(data$response$means))
    tran[[1+length(tran)]] = cbind(data$transition, UID = cur_uids[i])
    perf[[1+length(perf)]] = data.frame(nStates = data$nStates, MAPE = data$MAPE, R2 = data$R2, UID = cur_uids[i])    
    sdev[[1+length(sdev)]] = cbind(data.frame(t(data$response$stdev)), UID = cur_uids[i])
    names(sdev[[length(sdev)]]) = c(1:K, 'UID')
    
    # hack to account for those models with too few states
    if (is.list(data$response$stderr)){
      idx = which(as.numeric(lapply(data$response$stderr, function(x) which(is.na(x)))) == 1)
      err  = data$response$stderr
      err[[idx]] = data$response$stderr[[1]]
      err[[idx]] = err[[idx]] * 0
      err = as.data.frame(t(do.call('rbind', err)))        
    } else {
      err = data$response$stderr
    }

    if (ncol(err) < K) {
      missing_state = T
      idx = (ncol(err)+1):K
      err = cbind(err, rep(0, K - ncol(err)))  
    } else missing_state = F
    
    colnames(err) = 1:ncol(err)    
    serr[[1+length(serr)]] = cbind(as.data.frame(err), variable = rownames(err), UID = cur_uids[i])
    
    load(int_files[i])
    # hack to account for those models with too few states
    pi0 = data$benchmarks$pi0
    if (is.list(data$response$stderr) | missing_state) {
      pi0 = append(pi0, 0, after = idx-1)
      names(pi0) = 1:length(pi0)
    }
    init[[1+length(init)]] = pi0
  }
  init = as.data.frame(do.call('rbind', init))
  init = cbind(init, UID = cur_uids)  
  sdev = as.data.frame(do.call('rbind', sdev))
  serr = as.data.frame(do.call('rbind', serr))
  resp = as.data.frame(do.call('rbind', resp))
  tran = as.data.frame(do.call('rbind', tran))
  perf = as.data.frame(do.call('rbind', perf))
  
  return(list(init = init, sdev = sdev, serr = serr, resp = resp, tran = tran, perf = perf ))
}

# _____________________________________________________
# Constructor method for class DataImporter

setMethod(f = "initialize", 
          signature = "DataImporter",
          definition = function(.Object, 
                                path = './', 
                                verbose = T) {
            
            if (verbose) cat('*** Reading in data ***\n')
            
            # ____________________________
            # Format state attribute data
            
            res = read_rdata_files(path)

            .Object@RESPONSE    = res$resp
            .Object@RESP.ERR    = res$serr
            .Object@VOLATILITY  = res$sdev
            .Object@INIT.DIS    = res$init
            .Object@TRANSITION  = res$tran
            .Object@PERFORMANCE = res$perf
            
            return(.Object)
          })


# ____________________________________________________________
# Print method for class DataFormatter: display useful stats.

setMethod('show',
          signature  = 'DataImporter',
          definition = function(object){
            cat('*** DataImporter Object ***\n')
            
            # basic info
            cat(sprintf('No. Users:   %d\n', length(unique(object@INIT.DIS$UID))))            
            
            cat('*** END: DataImporter Object ***\n')
          })

# _________________________________________________________________
# Compute transition matrix given external covariates for one user

compute_transition_matrix = function(trans, x_var) {
  
  nStates= length(unique(trans$From))
  x_var  = cbind(rep(1,length(x_var)), x_var)
  prob   = lapply(1:nStates, function(s) {
    g = as.matrix(subset(trans, From == s))[,-which(colnames(trans) %in% c('From', 'To'))]
    y = x_var %*% t(g)
    y = exp(y) / rowSums(exp(y))
    return(y)
  })
  
  # build transition matrices
  prob = do.call('cbind', prob)
  dis = lapply(1:nrow(prob), function(i) {
    M = matrix(prob[i,], ncol = sqrt(ncol(prob)), byrow=T)
    return(M)              
  })

  return(dis)
}

# _____________________________________________________
# Propagate state probability distribution forward

propagate_distributions = function(p0, A_list){    
  P = list(as.matrix(p0) %*% A_list[[1]])  
  for (i in 2:length(A_list)) {
    A = A_list[[i]]
    P[[1+length(P)]] = P[[i-1]] %*% A
  }
  P = do.call('rbind', P)
  return(P)  
}

# ________________________________________________________
# Compute effective mean and variance of Gaussian mixture

compute_profile_statistics = function(p0, p, A_list, mu, sd){
  
  mu_eff = p %*% t(mu)
  sd_eff = sqrt(p %*% t(sd^2 + mu^2) - (p %*% t(mu))^2)
  covmat = matrix(0, nrow = nrow(p), ncol = nrow(p))
  for (t in 1:(nrow(covmat)-1)) {
    A = A_list[[t]]
    for (k in (t+1):ncol(covmat)) {
      A = A %*% A_list[[k]]
      covmat[t,k] = as.matrix(p0) %*% diag(mu) %*% A %*% t(mu) - (as.matrix(mu) %*% as.matrix(p[t,])) * (as.matrix(mu) %*% as.matrix(p[k,]))
    }
  }
  covmat = (covmat + t(covmat)) / 2
  diag(covmat) = sd_eff^2
  return(list(mu = mu_eff, sd = sd_eff, covmat = covmat))
}


# _________________________________________________
# Compute output profiles given a covariate series

setGeneric(
  name = "computeProfiles",
  def = function(.Object, UID_vec = NULL, x_vec = NULL){standardGeneric("computeProfiles")}
)
setMethod('computeProfiles',
          signature  = 'DataImporter',
          definition = function(.Object, UID_vec = NULL, x_vec = NULL) {
            
            # prepare inputs
            if (!is.null(UID_vec)) {
              UID_vec = intersect(.Object@INIT.DIS$UID, UID_vec)
            } else {
              UID_vec = as.character(.Object@INIT.DIS$UID)
            }                     
            if (is.null(x_vec)) {
              x_vec = rep(100, 24)
            }       
            
            # prepare distribution data
            trans  = subset(.Object@TRANSITION, UID %in% UID_vec)
            initd  = subset(.Object@INIT.DIS, UID %in% UID_vec)
            respm  = subset(.Object@RESPONSE, UID %in% UID_vec)
            respm  = respm[,-which(colnames(respm) == 'variable')]
            respe  = subset(.Object@RESP.ERR, UID %in% UID_vec)
            respe  = respe[,-which(colnames(respe) == 'variable')]
            volat  = subset(.Object@VOLATILITY, UID %in% UID_vec)
            
            # compute profile statistics
            distr  = mclapply(1:length(UID_vec), mc.cores =5, mc.preschedule =F, mc.silent = F, 
                            FUN = function(i) {
              
              uid = UID_vec[i]
              cat(paste(i, '/', length(UID_vec), '...', sep = ''))
              if (i %% 20 == 0) cat('\n')
              cur_trans = subset(trans, UID == uid)
              cur_initd = subset(initd, UID == uid)
              cur_respm = subset(respm, UID == uid)
              cur_respe = subset(respe, UID == uid)
              cur_volat = subset(volat, UID == uid)
              mu_a   = cur_respm[2,]  
              mu_b   = cur_respm[1,]  
              sd_a   = cur_respe[2,]  
              sd_b   = cur_respe[1,]                                        
              mu_a$UID = NULL
              mu_b$UID = NULL
              sd_a$UID = NULL
              sd_b$UID = NULL
              cur_trans$UID = NULL
              cur_initd$UID = NULL
              cur_volat$UID = NULL
              
              # propagate distributions
              A_list    = compute_transition_matrix(cur_trans, x_vec)
              cur_distr = propagate_distributions(cur_initd, A_list)
              
              # compute moments p0, p, A_list, mu, sd
              cur_profile_a = compute_profile_statistics(cur_initd, cur_distr, A_list, mu_a, sd_a)
              cur_profile_b = compute_profile_statistics(cur_initd, cur_distr, A_list, mu_b, sd_b)
              cur_profile_v = compute_profile_statistics(cur_initd, cur_distr, A_list, cur_volat * 0, cur_volat)
              
              return(list(pi = cur_distr, a = cur_profile_a, b = cur_profile_b, sigma = cur_profile_v$sd))
            })            
            names(distr) = UID_vec
            
            return(distr)     
          })      
          
# _________________________
# Method to plot schedules 

setMethod('plot',
          signature  = 'DataImporter',
          definition = function(x, type = 'default', inputs = NULL){            
            
            # population performance metrics
            if (type == 'default' | type == 'performance') {
              
              df = subset(x@PERFORMANCE, select = c('R2', 'MAPE'))
              df = melt(df)
              
              plt = ggplot(df) + 
                geom_density(aes(x = value, color = variable), size=2)
              plt = plt + scale_x_continuous(limits = c(0, 1))
              plt = plt + facet_wrap( ~ variable, nrow=1, scales = 'free') + 
                theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      axis.text.y      = element_text(size=18), 
                      axis.text.x      = element_text(size=18, angle = 330),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      plot.title       = element_text(size=20),   
                      legend.position  = c(0.8, 0.5),
                      axis.ticks = element_blank()) + 
                ylab('pdf') + xlab('') + ggtitle("Cross-Validation Performance")       
              
              return(plt)
            }
            
            # propagated profile structure for selected users: means + variance
            if (type == 'profiles') {
              
              if (is.null(inputs)) return('null selection!')
              profs.a.mu = lapply(inputs, function(y) y$a$mu)
              profs.a.mu = do.call('cbind', profs.a.mu)
              colnames(profs.a.mu) = names(inputs)
              profs.a.mu  = cbind(Hour = 1:nrow(profs.a.mu), as.data.frame(profs.a.mu))
              profs.a.mu  = melt(profs.a.mu, id.vars = c('Hour'))
              names(profs.a.mu) = c('Hour', 'User', 'mu')
              
              profs.a.se = lapply(inputs, function(y) y$a$sd)
              profs.a.se = do.call('cbind', profs.a.se)
              colnames(profs.a.se) = names(inputs)
              profs.a.se  = cbind(Hour = 1:nrow(profs.a.se), as.data.frame(profs.a.se))
              profs.a.se  = melt(profs.a.se, id.vars = c('Hour'))
              names(profs.a.se) = c('Hour', 'User', 'sd')              
              
              profs.b.mu = lapply(inputs, function(y) y$b$mu)
              profs.b.mu = do.call('cbind', profs.b.mu)
              colnames(profs.b.mu) = names(inputs)
              profs.b.mu  = cbind(Hour = 1:nrow(profs.b.mu), as.data.frame(profs.b.mu))
              profs.b.mu  = melt(profs.b.mu, id.vars = c('Hour'))
              names(profs.b.mu) = c('Hour', 'User', 'mu')              
              
              profs.b.se = lapply(inputs, function(y) y$b$sd)
              profs.b.se = do.call('cbind', profs.b.se)
              colnames(profs.b.se) = names(inputs)
              profs.b.se  = cbind(Hour = 1:nrow(profs.b.se), as.data.frame(profs.b.se))
              profs.b.se  = melt(profs.b.se, id.vars = c('Hour'))
              names(profs.b.se) = c('Hour', 'User', 'sd')              
              
              profs.a      = merge(profs.a.mu, profs.a.se, by = c("Hour", "User"))
              profs.b      = merge(profs.b.mu, profs.b.se, by = c("Hour", "User"))
              profs        = rbind(cbind(profs.a, Profile = 'Resp. [kWh/F]'), cbind(profs.b, Profile = 'Base [kWh]'))

              plt = ggplot(profs, aes(x = Hour, y = mu, color = Profile))
              plt = plt + geom_point(size = 2) + geom_line(size=1.5)
              plt = plt + geom_errorbar(aes(ymax = mu + sd, ymin = mu - sd))        
              #facet_grid(species~year, scales = "free_y", drop = FALSE)
              plt = plt + facet_grid( Profile ~ User, scales = 'free_y')
              plt = plt + theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      strip.text.y     = element_text(size=18), 
                      axis.text.y      = element_text(size=18), 
                      axis.text.x      = element_text(size=18),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      plot.title       = element_text(size=22),   
                      legend.position  = "none", 
                      axis.ticks = element_blank()) + 
                xlab('Hour of Day') + ggtitle("Selected Profiles")       
              return(plt)
            }
            
            # propagated profile structure for selected users: covariance matrix
            if (type == 'covmat') {

              # response matrices 
              a.mat = lapply(inputs, function(y) y$a$covmat)
              pla = ggplot(melt(a.mat), aes(x=X1, y=X2)) +
                facet_wrap(~ L1, nrow = 1) +
                geom_tile(aes(fill=value)) +
                coord_equal()
              pla = pla + scale_fill_gradient2(low = "blue", mid = "white", high = "red", space = 'rgb')
              pla = pla + theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      axis.text.y      = element_text(size=18), 
                      axis.text.x      = element_text(size=18),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      legend.title     = element_text(size=18),
                      legend.text     = element_text(size=18),
                      plot.title       = element_text(size=22),   
                      axis.ticks = element_blank()) + 
                xlab('Hour of Day') + ylab('Hour of Day') + ggtitle("Response Covariance")       
              
              # baseloads
              b.mat = lapply(inputs, function(y) y$b$covmat)
              plb = ggplot(melt(b.mat), aes(x=X1, y=X2)) +
                facet_wrap(~ L1, nrow = 1) +
                geom_tile(aes(fill=value)) +
                coord_equal()
              plb = plb + scale_fill_gradient2(low = "blue", mid = "white", high = "red", space = 'rgb')
              plb = plb + theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      axis.text.y      = element_text(size=18), 
                      axis.text.x      = element_text(size=18),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      legend.title     = element_text(size=18),
                      legend.text     = element_text(size=18),
                      plot.title       = element_text(size=22),   
                      axis.ticks = element_blank()) + 
                xlab('Hour of Day') + ylab('Hour of Day') + ggtitle("Baseload Covariance")       
              
              return(list(response = pla, baseload = plb))
            }
            
            # propagated profile structure for selected users: covariance matrix
            if (type == 'distributions') {
              
              mat = lapply(1:length(inputs), function(i) {                
                z = as.data.frame(inputs[[i]]$pi)
                names(z) = 1:ncol(z)
                z$User = names(inputs)[i]
                z$Hour = 1:nrow(z)
                z = melt(z, id.vars = c('User', 'Hour'))
                names(z) = c('User', 'Hour', 'Regime', 'Probability')
                return(z)
              })
              mat = do.call('rbind', mat)
              
              pla = ggplot(mat, aes(x=Hour, y=Probability, color = Regime))
              pla = pla + geom_point(size=3) + geom_line(size=1.5)
              pla = pla + facet_wrap(~ User, nrow = 1)
              pla = pla + theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      axis.text.y      = element_text(size=18), 
                      axis.text.x      = element_text(size=18),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      legend.title     = element_text(size=18),
                      legend.text     = element_text(size=18),
                      plot.title       = element_text(size=22),   
                      axis.ticks = element_blank())
              pla = pla + xlab('Hour of Day') + ylab('Pr[Regime(t)]') + ggtitle("Forecasted State Distributions")     
              
              return(pla)
            }
            
            # state space for selected users
            if (type == 'state-space') {
              
              # uid to name key
              tmp= names(inputs)
              names(tmp) = inputs
              
              # form plotting dataset
              df = subset(x@RESPONSE, UID %in% inputs)
              levels(df$variable) <- c('Baseload', 'Response')
              df$UID = tmp[as.character(df$UID)]
              df = melt(df)
              names(df)[3] = 'Regime'
              de = subset(x@RESP.ERR, UID %in% inputs)
              levels(de$variable) <- c('Base.Err', 'Resp.Err')
              de$UID = tmp[as.character(de$UID)]
              de = melt(de)
              names(de)[3] = 'Regime'
              dh = subset(x@VOLATILITY, UID %in% inputs)
              dh = melt(dh, id.vars = 'UID')
              names(dh)[2] = 'Regime'
              dh$variable = 'Volatility'
              dh = dh[,names(df)]
              dh$UID = tmp[as.character(dh$UID)]
              df = rbind(df, de, dh)
              df = cast(df, Regime + UID ~ variable)
              
              # compose plot
              pla = ggplot(df, aes(x = Baseload, y = Response, color = Regime))
              pla = pla + facet_wrap(~ UID, nrow = 1)
              pla = pla + geom_point(aes(size = Volatility)) + scale_size_area(max_size = 12) 
              pla = pla + geom_errorbar(aes(ymax = Response + Resp.Err, ymin = Response - Resp.Err))                            
              pla = pla + geom_errorbarh(aes(xmax = Baseload + Base.Err, xmin = Baseload - Base.Err))                            
              pla = pla + theme_bw() + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      strip.text.x     = element_text(size=18),
                      axis.text.y      = element_text(size=18), 
                      axis.text.x      = element_text(size=18),
                      axis.title.y     = element_text(size=18),
                      axis.title.x     = element_text(size=18),
                      legend.title     = element_text(size=18),
                      legend.text     = element_text(size=18),
                      plot.title       = element_text(size=22),   
                      axis.ticks = element_blank()) + 
                xlab('Baseload [kWh]') + ylab('Response [kWh/F]') + ggtitle("State Space Characterization")       
              
              return(pla)
            }
          })


