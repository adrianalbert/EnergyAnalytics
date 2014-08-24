# StateDecoder.r
#
# Extracts occupancy states from obs time series data.
# 
# Adrian Albert
# Last modified: February 2014.
# -----------------------------------------------------------------------

library('timeDate')
require('lmtest')
library('depmixS4')
library('R.utils')
library('MASS')
library('lubridate')
library('dummies')

source('./estimation/viterbi_states.R')

### Remove !!! 
# Hack into the depmix direct likelihood maximization
# source('./estimation/depmixfit-modified.R')
# source('classes/responseTruncNORM.R')
### Remove !!! 

# clean-up previous definitions of methods for class StateDecoder
removeClass('StateDecoder')

options(error = recover)

# ________________________
# Class definition
setClass(
  Class = "StateDecoder",
  representation = representation(
    UID        = "character",         # unique person ID
    data.train = 'data.frame',
    data.test  = 'data.frame',
    resp.vars  = 'character',         # covariates for analysis (response)
    tran.vars  = 'character',         # covariates for analysis (transition)
    controls   = 'list',              # controls for fitting the model
    HMM        = "list",               # summary of HMM model
    performance= "list"
    )
)

# _______________________________________
# Constructor method for class StateDecoder

setMethod(f = "initialize", 
          signature = "StateDecoder",
          definition = function(.Object, data, timestamps, UID, 
                                train.frac = 1, 
                                tran.vars = c(), resp.vars = c(), 
                                controls = list(
                                  Kmin = 2, Kmax = 4, 
                                  maxit = 100, nRestarts = 5,
                                  thresh.R2 = 0.85, thresh.MAPE = 0.15), 
                                verbose = T) {
            
            .Object@UID = UID
            
            if (verbose) {
              cat(paste('*** Initializing StateDecoder (', .Object@UID, ') ***\n', sep=''))
              t0 = proc.time()
            }
                        
            # retain only interesting data
            data = subset(data, select = c('obs', setdiff(resp.vars, '(Intercept)'), setdiff(tran.vars, '(Intercept)')))
            
            # add in timestamps back
            data$timestamps = timestamps
            
            # split up into test and train
            T_train   = trunc(nrow(data) * train.frac)
            # T_train   = trunc(T_train / 24) * 24
            idx_train = 1:T_train
            idx_test  = setdiff(1:nrow(data), idx_train)    
            if (length(idx_train) > 0) data.train = data[idx_train,] else data.train = data.frame()                                      
            if (length(idx_test) > 0) data.test = data[idx_test,] else data.test = data.frame()                                      
            
            # set response and transition covariates
            if (!is.null(resp.vars)) .Object@resp.vars = resp.vars
            if (!is.null(tran.vars)) .Object@tran.vars = tran.vars
            
            .Object@data.train = data.train
            .Object@data.test  = data.test
            .Object@controls   = controls
            
            .Object@HMM = list()
            
            if (verbose){
              dt = proc.time() - t0
              print(dt)
            }
            
            return(.Object)
          })

# ____________________________________________________________
# Print method for class StateDecoder: display useful stats.

setMethod('show',
          signature  = 'StateDecoder',
          definition = function(object){
            cat('*** StateDecoder Object ***\n')
            
            # basic info
            cat(sprintf('UID:         %s\n', object@UID))            
            if (length(object@HMM) > 0) {
      	      cat(sprintf('HMM States:\t%d\n', object@HMM$nStates))
      	      cat(sprintf('HMM MAPE =\t%f; R2 =\t%f\n', object@HMM$MAPE, object@HMM$R2))
      	      cat("Response\n")
      	      print(object@HMM$response)
      	      cat("Transition\n")
      	      print(object@HMM$transition)
            }
            
            cat('*** END: StateDecoder Object ***\n')
          })

# ______________________________________
# Define HMM model for depmixS4 package

defineHMM.depmix = function(data.train, K = 3, type = 'default', 
                               response.vars = NULL, transitn.vars = NULL, 
                               respstart = NULL, trstart = NULL) {
  
  intercept_resp = '(Intercept)' %in% response.vars
  intercept_tran = '(Intercept)' %in% transitn.vars
  response.vars  = setdiff(response.vars, '(Intercept)')
  transitn.vars  = setdiff(transitn.vars, '(Intercept)')
  
  # define response model
  if (intercept_resp) fmla_response = 'obs ~ 1' else fmla_response = 'obs ~ -1'
  if (length(response.vars)>0) fmla_response = paste(fmla_response, paste(response.vars, collapse = '+'), sep='+')
  fmla_response = as.formula(fmla_response)            
  
  # define transition model
  if (intercept_tran) fmla_transitn = ' ~ 1' else fmla_transitn = ' ~ -1'
  if (length(transitn.vars)>0) fmla_transitn = paste(fmla_transitn, paste(transitn.vars, collapse = '+'), sep='+')
  fmla_transitn = as.formula(fmla_transitn)   
      
  # set up model for depmixS4
  if (type == 'default') {
      
#     mod <- makeModelTruncNORM(data.train, fmla_response, fmla_transitn, K = K)
      
      mod <- depmix(response  = fmla_response, 
                  transition= fmla_transitn,
                  data      = data.train, 
                  nstates   = K, 
                  respstart = respstart,
                  trstart   = trstart,
                  family    = gaussian(),
                  instart   = runif(K))
  } else {
    mod = makeModelSelect(data.train, fmla_response, fmla_transitn, K = K)    
  }
  
  return(mod)
}



# ______________________________________
# Define constraints on the model

defineConstraints = function(mod) {
  
  pars  = c(unlist(getpars(mod)))
  fixed = rep(0, length(pars))
  K_opt = nstates(mod) 
  nResp = length(mod@response[[1]][[1]]@parameters$coefficients)
  
	# define non-negative base levels
  conMat= matrix(0, ncol = npar(mod), nrow=K_opt)
  bl    = rep(0, K_opt)
  bu    = rep(20, K_opt)
  for (j in 1:K_opt) {
    idx.base = K_opt + K_opt^2 * nResp + (nResp + 1) * (j-1) + 1
    conMat[j,idx.base] = 1
  }
  
  # define 4-state model
	# state 1 is bursty (B)
	for (j in 1:K_opt){
	  pars[K_opt*j*nResp + 1] = 0 # no temperature dependence in transition
	}
  pars[K_opt + K_opt^2*nResp + nResp] = 0 # no temperature dependence in response
	fixed[K_opt*1:K_opt*nResp + 1] = 1
	fixed[K_opt + K_opt^2*nResp + nResp] = 1
	
	# state 2 is no-HVAC (N)
  pars[K_opt + K_opt^2*nResp + nResp*2 + 1] = 0 # no temperature dependence in response
	fixed[K_opt + K_opt^2*nResp + nResp*2 + 1] = 1
	
	mod.cs = setpars(mod, pars)
  return(list(conMat = conMat, lower = bl, upper = bu, mod = mod.cs, fixed = fixed))  
}

# ______________________________________
# Wrapper to fit model

fit.model = function(mod, maxit = 100, tol = 1e-3){ 
  
   # constr = defineConstraints(mod)
    mod.cs = mod #constr$mod

		# for donlp
    # contrl = donlp2Control()
    # contrl$epsx = 1e-4
    # contrl$epsfcn = 1e-8
    # contrl$silent = F
    # contrl$te1 <- contrl$te2 <- contrl$te3 <- T  
    # contrl$nreset.multiplier = 2

    # print(constr)
    
    print(tol)
  
    mod.fit = fit(mod.cs, verbose = T, emcontrol = em.control(maxit = maxit, tol = tol))
    
#    mod.fit = fit(mod.cs, verbose = T, 
#                  conrows = constr$conMat, 
#                  conrows.lower = constr$lower, 
#                  conrows.upper = constr$upper, 
#                  fixed  = constr$fixed,
#                  method = 'rsolnp',
#                  solnpcntrl=list(rho = 1, outer.iter = 100, inner.iter = 200, delta = 1e-5, tol = 1e-6, trace=1))
#                  donlpcntrl=contrl)

  return(mod.fit)
}

# ______________________________________
# Wrapper that treats for fitting errors

fitHMM = function(mod, nRestarts = 1, verbose = T, maxit = 100, tol = 1e-3){ 

  if (verbose) cat('---> Learning HMM...\n')
    
  # fit given model            
  ok = FALSE
  it = 0
  cur_tol = tol
  
  while (!ok & it <= nRestarts) {
    
#     out <- capture.output(fm  <- try(fit.model(mod, maxit = maxit, cur_tol = cur_tol)))                
#     nlines = length(out)
    fm  <- try(fit.model(mod, maxit = maxit, tol = cur_tol))
    nlines = 3
    
    if (class(fm) != 'try-error') {
      ok = nlines > 2                                   
    } else ok = FALSE
    if (!ok){
      if (verbose) cat('Bad HMM fit; re-estimating...\n')
      cur_tol = cur_tol / 10
    } else {
      if (verbose) cat(paste('Convergence in', nlines*5, 'iterations.\n'))
    }
    it = it + 1
  }
  return(fm)  
}            

# __________________________________________
# Cross-validation for HMM fitting.

fitHMM.cv = function(data, K = 3, 
                     response.vars = NULL, transitn.vars = NULL, 
                     respstart = NULL, trstart = NULL, 
                     nRestarts = 1, verbose = T, maxit = 100, tol = 1e-3){ 
  if (verbose)
    cat(paste('---> HMM Cross-Validation K =', K,'\n'))
  
  # fit HMM on 1/2 of data
  data_1  = data[seq(1, nrow(data)-1, by=2),]
  mod     = defineHMM.depmix(data_1, K = K, 
                      response.vars = response.vars, 
                      transitn.vars = transitn.vars, 
                      respstart     = respstart,
                      trstart       = trstart)

  fm_half = fitHMM(mod, nRestarts = nRestarts, verbose = verbose, maxit = maxit, tol = tol)
  
  if (class(fm_half) == 'try-error') 
    return(list(model.half = NA, metrics = c(MAPE.cv = NA, R2.cv = NA, BIC = NA, AIC = NA)))
  
  # decode the other half
  data_2  = data[seq(2, nrow(data), by=2),]
  data.dum= dummy.data.frame(data_2[,-which(names(data_2)=='timestamps')])
  probs   = viterbi_states(fm_half, data.dum, data.dum$obs)
  
  # compute Cross-Validation error metrics
  MAPE    = mean(abs((probs$fit - data_2$obs) / probs$fit))
  R2      = 1 - sum((probs$fit - data_2$obs)^2) / sum((data_2$obs - mean(data_2$obs))^2)
  
  return(list(model.half = fm_half, metrics = c(MAPE.cv = MAPE, R2.cv = R2, BIC = BIC(fm_half), AIC = AIC(fm_half))))
}            

# _________________________________
# Learn HMM and optimum model size

learnModelSize = function(data.train, resp.vars = c(), tran.vars = c(), 
                             verbose = T, 
                             Kmin = 3, Kmax = 3, 
                             nRestarts = 1, maxit = 100, tol = 1e-3, 
                             thresh.R2 = 0.8, thresh.MAPE = 0.15) { 

    # _______________________________________________
    # choose model size K (number of states)

    if (Kmax > Kmin) {
      result = list()
      k      = Kmin
      done   = F
      while (k <= Kmax & !done) {
        result[[as.character(k)]] = fitHMM.cv(data.train, K = k, 
                                              response.vars = resp.vars,
                                              transitn.vars = tran.vars, 
                                              respstart = NULL, trstart = NULL, 
                                              nRestarts = nRestarts, verbose = verbose, maxit = maxit, tol = tol)
        R2.cv   = result[[as.character(k)]]$metrics['R2.cv']
        MAPE.cv = result[[as.character(k)]]$metrics['MAPE.cv']
        if (R2.cv > thresh.R2 | MAPE.cv < thresh.MAPE) done = T else k = k + 1
      }
      metric = t(sapply(result, function(l) l$metrics))
      rownames(metric) = names(result)
      if (done) K_opt = k else {            
        # get minimum BIC from good fits 
        BICs   = sapply(result, function(l) l$metrics[3])
        bad.fit= which(!is.finite(BICs))
        if (length(BICs) == length(bad.fit)) {
          stop('HMM fitting error!')
        }
        Kvec   = Kmin:Kmax
        idx_opt= which.max(metric[,'R2.cv'])             
        K_opt  = Kvec[idx_opt]
      }
    } else { 
      K_opt = Kmin   
      metric = NULL
      result = NULL
    }
                          
    # _______________________________________________
    # Fit model to full data
    
    trstart = NULL

    # fit model to full data 
    mod    = defineHMM.depmix(data.train, K = K_opt, 
                       response.vars = resp.vars, 
                       transitn.vars = tran.vars, 
                       respstart = NULL, 
                       trstart = trstart)
    fm     = fitHMM(mod, nRestarts = nRestarts, verbose = verbose, maxit = maxit, tol = tol)
    
    # ______________________
    # Save computation
    
    rm(list = c('result'))
    gc()
    return(list(model = fm, size = K_opt, metrics = metric))
}           
         
# _______________________________________________
# Compute estimates of out-of-sample performance

computePredictionAccuracy = function(model, data.test, test.periods = 5) {
            
    TT  = nrow(data.test)
    per = (1:TT) %/% (TT %/% test.periods) + 1
    if (length(which(per == test.periods+1))<=2) per[which(per == test.periods+1)] = test.periods
    test.periods = length(unique(per))
    
    data.dum = dummy.data.frame(data.test[,-which(names(data.test)=='timestamps')])
    
    stats = list()
    for (p in 1:test.periods) {

      idx_test   = which(per <= p) 
      T_test     = length(idx_test)
      test_data  = data.dum[idx_test,]                        
      
      if (verbose) cat(paste('Test window:', T_test, '\n'))

      # set up Viterbi decoding on test data
      probs      = viterbi_states(model, test_data, test_data$obs)              
      
      # compute decoding performance statistics
      MAPE       = mean(abs((probs$fit - test_data$obs) / probs$fit))
      R2         = 1 - sum((probs$fit - test_data$obs)^2) / sum((test_data$obs - mean(test_data$obs))^2)
                    
      # store stats
      stats[[as.character(T_test)]] = c(MAPE, R2)
    }
    len_vec = names(stats)
    stats   = do.call('rbind', stats)
    rownames(stats) = len_vec
    colnames(stats) = c('MAPE', 'R2')
    
    return(stats)
}

# ________________________________________
# Compute statistics for HMM fit

# function to compute Generalized R2
# GR2 = function(LL0, LL, N) {
#   gr2 = 1 - exp(2*(LL0 - LL) / N)
# }

# function to compute R2
R2 = function(y, yfit) {
  r2 = 1 - sum((y - yfit)^2) / sum((y-mean(y))^2)
}

# ________________________________________________
# Extract useful parameters out of depmixS4 model

extractParameters.depmixS4 = function(fm_opt){
              
  HMM   = list()
  K_opt = length(fm_opt@response)
  HMM[['nStates']] = K_opt
  
  # ________________________________
  # Response parameters
  
  resp_vars  = colnames(fm_opt@response[[1]][[1]]@x)
  HMM[['response']] = list()
  stddev     = sapply(1:K_opt, function(j) fm_opt@response[[j]][[1]]@parameters$sd)
  names(stddev) = 1:K_opt
  HMM[['response']][['stdev']] = stddev
  params = sapply(1:K_opt, function(j) fm_opt@response[[j]][[1]]@parameters$coefficients)
  if (class(params) == 'numeric') params = data.frame(params) else params = data.frame(t(params))
  rownames(params) = 1:K_opt
  colnames(params) = resp_vars
  HMM[['response']][['means']]  = as.data.frame(t(params))
  
  # ______________________________
  # Transition parameters

  # format into matrix for output 
  transitn_vars  = colnames(fm_opt@transition[[1]]@x)          
  params = lapply(1:K_opt, function(j) {
    tmp = fm_opt@transition[[j]]@parameters$coefficients
    rownames(tmp) = transitn_vars              
    if (class(tmp) != 'numeric') return(data.frame(t(tmp))) else return(data.frame(tmp))
  })                          
  params = do.call(rbind, lapply(params, data.frame, stringsAsFactors=FALSE))     
  params$From = rep(1:K_opt, each = K_opt)
  params$To   = rep(1:K_opt, times = K_opt)
  if (!is.null(ncol(params))) colnames(params) = c(transitn_vars, 'From', 'To')
  HMM[['transition']] = params
  
  # Viterbi states
  HMM[['states']] = posterior(fm_opt)
  
  # ________________________________________
  # Perform preliminary analysis of errors
  
  # Predicted state means
  fit = sapply(1:K_opt, function(k) predict(fm_opt@response[[k]][[1]])) 
  HMM[['fit']] = sapply(1:nrow(fit), function(j) fit[j,HMM[['states']][j,1]])
  HMM[['fit.avg']]     = sapply(1:nrow(fit), function(j) sum(fit[j,] * HMM[['states']][j,-1]))
    
  # model fit (penalized likelihood)
  HMM[['BIC']]    = BIC(fm_opt)  
              
  rm(list = c('fm_opt'))
  gc()
  return(HMM)
}

# ________________________________________
# Compute decoding performance statistics 

estimateStandardErrors = function(states, data, resp.vars) {
  
  K = length(unique(states))
  resp.vars = setdiff(resp.vars, '(Intercept)')
  stderr = sapply(1:K, function(k) {
    fmla = as.formula(paste('obs ~ ', paste(resp.vars, collapse = '+')))
    fit  = lm(fmla, data[which(states == k),])
    se   = summary(fit)$coefficients[,2]
    return(se)
  })
  
  return(stderr)
  
}

# ________________________________________
# Compute decoding performance statistics 

setGeneric(
  name = "computeDecodingStats",
  def = function(.Object, verbose = T){standardGeneric("computeDecodingStats")}
)
setMethod('computeDecodingStats',
          signature  = 'StateDecoder',
          definition = function(.Object, verbose = T){
            
            # residuals for each state
            residuals = .Object@data.train$obs - .Object@HMM$fit                          
            .Object@HMM[['residual']] = residuals
            
            # test normality of residuals in each state
            is_normal = sapply(1:.Object@HMM[['nStates']], function(j) {
              idx = which(.Object@HMM[['states']][,1] == j)
              if (length(na.omit(residuals[idx]))>3 & length(na.omit(residuals[idx]))<5000)
                res = shapiro.test(na.omit(residuals[idx])) else res = list(p.value=NA)
              pval= res$p.value
              return(pval)
            })
            .Object@HMM$sw_test = is_normal
            
            # in-sample fit performance metrics
            .Object@HMM[['MAPE']]    = mean(abs(residuals / .Object@HMM$fit))
            .Object@HMM[['R2']]      = R2(.Object@data.train$obs, .Object@HMM$fit)
            
            return(.Object)
          })

# _______________________________________________
# Populate StateDecoder object

setGeneric(
  name = "learnStateDecoder",
  def = function(.Object, verbose = T){standardGeneric("learnStateDecoder")}
)
setMethod('learnStateDecoder',
          signature  = 'StateDecoder',
          definition = function(.Object, verbose = T){
            
            if (verbose) {
              cat(paste('*** HMM analysis for StateDecoder UID: (', .Object@UID, ')***\n', sep=''))
              t0 = proc.time()
            }
            
            # learn model 
            controls = .Object@controls
            model = learnModelSize(.Object@data.train, 
                                   resp.vars = .Object@resp.vars, tran.vars = .Object@tran.vars, 
                                   verbose = T, 
                                   Kmin = controls$Kmin, Kmax = controls$Kmax, 
                                   nRestarts = controls$nRestarts, maxit = controls$maxit, tol = controls$tol, 
                                   thresh.R2 = controls$thresh.R2, thresh.MAPE = controls$thresh.MAPE)
            
            # compute prediction accuracy out-of-sample
            performance          = list()            
            performance$accuracy = NULL #computePredictionAccuracy(model$model, .Object@data.test, test.periods = 5)
            performance$cv.stats = model$metrics
            .Object@performance  = performance
                            
            # compute decoding performance stats
            .Object@HMM = extractParameters.depmixS4(model$model)
                                        
            # estimate state-specific coefficient standard errors
            # note that this over-estimates errors because all variance is attributed just to errors in state-based
            # coefficients, not also to the coefficients of the transition matrix
            .Object@HMM$response$stderr = estimateStandardErrors(.Object@HMM$states[,1], .Object@data.train, .Object@resp.vars)
                        
            # compute some stats
            .Object = computeDecodingStats(.Object)
            
            if (verbose) {
              dt = proc.time() - t0;
              print(dt)
            }
            
            return(.Object)
              
          })

# ______________________________________
# Method to dump StateDecoder model data

setGeneric(
  name = "dumpDecodedData",
  def = function(.Object, verbose = T, path = NULL){standardGeneric("dumpDecodedData")}
)
setMethod('dumpDecodedData',
          signature  = 'StateDecoder',
          definition = function(.Object, verbose = T, path = NULL){
            if (verbose) 
              cat(paste('*** Dumping decoded data (', .Object@UID, ')***\n', sep=''))
            
            data          = .Object@HMM
            data$UID      = .Object@UID
            data$states   = data$states[,1]
            data$fit      = NULL
            data$residual = NULL
            data$fit.avg  = NULL
            
            if (is.null(path)) return(data) else {
              save(list = c('data'), file = paste(path, .Object@UID, '_decoded.RData', sep=''))
              return(NULL)
            }
          })
