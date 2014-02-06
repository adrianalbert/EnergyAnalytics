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

# clean-up previous definitions of methods for class 
removeClass('DataImporter')

# ________________________
# Class definition

setClass(
  Class = "DataImporter",
  representation = representation(
    RESPONSE        = "data.frame",          # state parameters
    RESP.ERR        = "data.frame",          # state parameters errors
    INIT.DIS        = "data.frame",          # initial distribution over states
    TRANSITION      = "data.frame",          # transition parameters
    PERFORMANCE     = "data.frame"           # transition parameters
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
  sder       = list()
  init       = list()
  tran       = list()
  perf       = list()  
  for (i in 1:length(cur_uids)) {
    load(dec_files[i])
    resp[[1+length(resp)]] = cbind(data$response$means, UID = cur_uids[i], variable = rownames(data$response$means))
    tran[[1+length(tran)]] = cbind(data$transition, UID = cur_uids[i])
    perf[[1+length(perf)]] = data.frame(nStates = data$nStates, MAPE = data$MAPE, R2 = data$R2, UID = cur_uids[i])    
    K = ncol(data$response$means)
    
    # hack to account for those models with too few states
    if (is.list(data$response$stderr)){
      idx = which(as.numeric(lapply(data$response$stderr, function(x) which(is.na(x)))) == 1)
      err  = data$response$stderr
      err[[idx]] = data$response$stderr[[1]]
      err[[idx]] = err[[idx]] * 0
      err = as.data.frame(t(do.call('rbind', err)))  
      colnames(err) = 1:ncol(err)
      err = err[2,]
    } else err = data$response$stderr[2,]
    if (length(err) < K) {
      missing_state = T
      idx = (length(err)+1):K
      err = c(err, rep(0, K - length(err)))  
    } else missing_state = F
    sder[[1+length(sder)]] = cbind(as.data.frame(rbind(stdev = data$response$stdev, stderr = err)), 
                                   variable = c('stdev', 'sderr'),
                                   UID = cur_uids[i])
    
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
  sder = as.data.frame(do.call('rbind', sder))
  resp = as.data.frame(do.call('rbind', resp))
  tran = as.data.frame(do.call('rbind', tran))
  perf = as.data.frame(do.call('rbind', perf))
  
  return(list(init = init, sder = sder, resp = resp, tran = tran, perf = perf ))
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
            .Object@RESP.ERR    = res$sder
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
  sd_eff = p %*% t(sd^2 + mu^2) - (p %*% t(mu))^2
  covmat = matrix(0, nrow = nrow(p), ncol = nrow(p))
  for (t in 1:(nrow(covmat)-1)) {
    A = A_list[[t]]
    for (k in (t+1):ncol(covmat)) {
      A = A %*% A_list[[k]]
      covmat[t,k] = as.matrix(p0) %*% diag(mu) %*% A %*% t(mu) - (as.matrix(mu) %*% as.matrix(p[t,])) * (as.matrix(mu) %*% as.matrix(p[k,]))
    }
  }
  covmat = (covmat + t(covmat)) / 2  
  diag(covmat) = sd_eff
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
            
            # compute profile statistics
            distr  = mclapply(1:length(UID_vec), mc.cores =5, mc.preschedule =F, mc.silent = F, FUN = function(i) {
              
              uid = UID_vec[i]
              cat(paste(i, '/', length(UID_vec), '...', sep = ''))
              if (i %% 20 == 0) cat('\n')
              cur_trans = subset(trans, UID == uid)
              cur_initd = subset(initd, UID == uid)
              cur_respm = subset(respm, UID == uid)
              cur_respe = subset(respe, UID == uid)
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
              
              # propagate distributions
              A_list    = compute_transition_matrix(cur_trans, x_vec)
              cur_distr = propagate_distributions(cur_initd, A_list)
              
              # compute moments p0, p, A_list, mu, sd
              cur_profile_a = compute_profile_statistics(cur_initd, cur_distr, A_list, mu_a, sd_a)
              cur_profile_b = compute_profile_statistics(cur_initd, cur_distr, A_list, mu_b, sd_b)
              
              return(list(pi = cur_distr, a = cur_profile_a, b = cur_profile_b))
            })            
            names(distr) = UID_vec
            
            return(distr)     
          })      
          
# _________________________
# Method to plot schedules 

setMethod('plot',
          signature  = 'DataImporter',
          definition = function(x, type = 'default', selected = NULL){            
            
            # population performance metrics
            if (type == 'default' | type == 'performance') {
              
              df.1 = subset(x@OCCUP_STATS, select = c('R2.cv', 'MAPE.cv'))
              names(df.1) = c('R2', 'MAPE')
              df.1 = melt(df.1)
              df.1$Type = 'Out.Of.Sample'
              df.2 = subset(x@OCCUP_STATS, select = c('R2', 'MAPE'))
              df.2 = melt(df.2)
              df.2$Type = 'In.Sample'
              df = rbind(df.1, df.2)
              
              plt = ggplot(df) + 
                geom_density(aes(x = value, colour = Type), size=2)
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
                ylab('pdf') + ggtitle("Occupancy States: Cross-Validation Performance")       
              
              return(plt)
            }
            
            # propagated profile structure for selected users: means + variance
            
            # propagated profile structure for selected users: covariance matrix
            
            
            
          })


