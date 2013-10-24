# responseBreakpointFree.r
#
# Define breakpoint model for energy with temperature.
# 
# Adrian Albert
# Last modified: April 2013.
# -----------------------------------------------------------------------

library('segmented')

# define a response class which only contains the standard slots, no additional slots
setClass("breakpointFree", 
         representation(formula="formula",
                        family="ANY",
                        state_type = 'character'),
         contains="response")

# define a generic for the method defining the response class
setGeneric("breakpointFree",
           function(formula, breakpoint.var = NULL, data = NULL, pstart = NULL, fixed = NULL, prob=TRUE, na.action="na.pass", ...) standardGeneric("breakpointFree"))

# define the method that creates the response class
setMethod("breakpointFree", 
          signature(formula = 'formula'), 
          function(formula, breakpoint.var = NULL, data=NULL, pstart=NULL, fixed=NULL, prob=TRUE, na.action="na.pass", ...) {
            
            # parse input parameters
            call <- match.call()
            mf <- match.call(expand.dots = FALSE)
            m <- match(c("formula", "data"), names(mf), 0)
            mf <- mf[c(1, m)]
            mf$drop.unused.levels <- TRUE
            mf$na.action <- na.action
            mf[[1]] <- as.name("model.frame")
            mf <- eval(mf, parent.frame())
            x <- model.matrix(attr(mf, "terms"),mf)
            if(any(is.na(x))) stop("'depmixS4' does not currently handle covariates with missing data.")
            y <- model.response(mf)
            if(!is.matrix(y)) y <- matrix(y,ncol=1)
            colnames(y) = all.vars(formula)[1]
                        
            # if no breakpoint variable is defined assume first variable in formula
            if (is.null(breakpoint.var)) breakpoint.var = all.vars(formula)[2]
            
            coefs.init = vector("numeric",length=ncol(x))
            
            # set output parameters
            parameters <- list()
            constr <- NULL
            parameters$coefficients <- coefs.init
            parameters$sd <- 1
            parameters$psi = mean(x[,breakpoint.var])
            parameters[[breakpoint.var]] = 0            
            constr <- list(
              parup = rep(Inf,length(parameters$coefficients)+3),
              parlow = c(rep(-Inf,length(parameters$coefficients)+2),0)
            )
            
            npar <- length(unlist(parameters))
            if(is.null(fixed)) fixed <- as.logical(rep(0,npar))
            if(!is.null(pstart)) {
              if(length(pstart)!=npar) stop("length of 'pstart' must be ",npar)
              parameters$coefficients[1:length(parameters$coefficients)] <- pstart[1:length(parameters$coefficients)]
              if(length(unlist(parameters))>length(parameters$coefficients)) {
                parameters$sd   <- as.numeric(pstart[(length(parameters$coefficients)+1)])
                parameters$psi  <- as.numeric(pstart[(length(parameters$coefficients)+2)])
                names(parameters$psi) = breakpoint.var
                parameters[[breakpoint.var]] = 0 #as.numeric(pstart[(length(parameters$coefficients)+3)])                
              }
            }                    
            
            # mod <- new("exgaus",parameters=parameters,fixed=fixed,x=x,y=y,npar=npar)
            mod = new(Class = "breakpointFree", formula = formula, state_type = '',
                      parameters=parameters, fixed=fixed, x=x,y=y,npar=npar,constr=constr)
            # mod@state_type = ''
          }
)

setMethod("show","breakpointFree",
          function(object) {
            cat("Model of type breakpoint (see ?segmented for details) formula: ", sep="")
            print(object@formula)            
            cat(sprintf("Type: %s\n",object@state_type))
            cat("Coefficients: \n")
            print(object@parameters$coefficients)
            cat("sd ",object@parameters$sd,"\n")
            cat(paste('Break point', names(object@parameters$psi),'\n'))
            print(object@parameters$psi)
            cat('Difference in slope\n')
            print(object@parameters[[names(object@parameters$psi)]])            
          }
)

setMethod("predict","breakpointFree",
          function(object) {
            break.var = names(object@parameters$psi)
            idx       = which(colnames(object@x) == break.var)
            if (ncol(object@x) >= 2) {
              y.pred = object@x %*% object@parameters$coefficients
            } else 
              y.pred = object@x * object@parameters$coefficients
            indict = (object@x[,idx] - object@parameters$psi) > 0
            y.adj  = (object@x[,idx] - object@parameters$psi) * object@parameters[[break.var]] * indict
            y.pred = y.pred + y.adj
            return(y.pred)
          }
)

setMethod("dens","breakpointFree",
          function(object,log=FALSE) {
            dnorm(x=object@y,mean=predict(object),sd=object@parameters$sd,log=log)
          }
)

setGeneric("logDens",function(object,...) standardGeneric("logDens"))
setMethod("logDens","breakpointFree",
          function(object) {
            dnorm(x=object@y,mean=predict(object),sd=object@parameters$sd,log=TRUE)
          }
)

setMethod("logLik","breakpointFree",
          function(object) {
            sum(logDens(object))
          }
)

setMethod("simulate",signature(object="breakpointFree"),
          function(object,nsim=1,seed=NULL,times) {
            if(!is.null(seed)) set.seed(seed)
            if(missing(times)) {
              # draw in one go
              mu <- predict(object)
            } else {
              mu <- predict(object)[times]
            }  
            nt <- length(mu)
            sd <- object@parameters$sd
            response <- rnorm(nt*nsim,mean=mu,sd=sd)
            #if(nsim > 1) response <- matrix(response,ncol=nsim)
            response <- as.matrix(response)
            return(response)
          }
)

library('segmented')
setMethod("fit","breakpointFree",
          function(object, w, maxit = 50, tol = 1e-2) {
            
            if(missing(w)) w <- NULL
            nas <- is.na(rowSums(object@y))
            pars <- object@parameters
            break.var = names(object@parameters$psi)
            
            # first prepare data 
            data          = as.data.frame(cbind(object@x[!nas,], object@y[!nas,]))
            names(data)   = c(colnames(object@x), colnames(object@y))            
            idx_intercept = which(names(data) == '(Intercept)')
            if (length(idx_intercept)>0) data = data[,-idx_intercept]
            predictors    = intersect(names(data), colnames(object@x))

            # define models to be estimated
            if (length(predictors)>1) {
              fmla_1      = as.formula(paste(colnames(object@y), '~', paste(setdiff(predictors, break.var), collapse = '+')))
              fmla_2      = as.formula(paste(colnames(object@y), '~', paste(predictors, collapse = '+')))
            } else {
              fmla_1      = as.formula(paste(colnames(object@y), '~ 1'))
              fmla_2      = as.formula(paste(colnames(object@y), '~ ', break.var))
            }
            fmla     = as.formula(paste('~',break.var))
            
            # estimate mean-only model
            if(!is.null(w)) {              
              fit.1 <- lm(fmla_1, data = data, weights = w[!nas])
            } else {
              fit.1 <- lm(fmla_1, data = data)
            }
            
            # estimate slope model
            if(!is.null(w)) {              
              fit.2 <- lm(fmla_2, data = data, weights = w[!nas])
            } else {
              fit.2 <- lm(fmla_2, data = data)
            }

#             # estimate right-slope model (null left slope)
#             fit.seg.right  <- try(segmented(fit.1, seg.Z = fmla, psi = pars$psi))
#             
#             # estimate left-slope model (right slope null)
#             data.neg = data
#             data.neg[,break.var] = -data.neg[,break.var]
#             if(!is.null(w)) {              
#               fit.neg <- lm(fmla_1, data = data.neg, weights = w[!nas])
#             } else {
#               fit.neg <- lm(fmla_1, data = data.neg)
#             }            
#             fit.seg.left <- try(segmented(fit.neg, seg.Z = fmla, psi = -pars$psi))            
            
            # estimate dual-slope model
            fit.seg.dual <- try(segmented(fit.2, seg.Z = fmla, psi = pars$psi))
            
            # compute best model
            models = list(flat = fit.1, slope = fit.2, 
                          # left = fit.seg.left, right = fit.seg.right)
                           dual = fit.seg.dual)
            lglik_vec = sapply(models, function(fit) {
              if (class(fit)[1] == 'try-error') return(NA)
              return(BIC(fit))
            })
                        
            # compute variance for each model
            err_vec = sapply(models, function(fit) {
              if (class(fit)[1] == 'try-error') return(NA)
              if(!is.null(w)) {
                cur_sd <- sqrt(sum(w[!nas]*fit$residuals^2/sum(w[!nas])))
              } else {
                cur_sd <- sd(fit$residuals)
              }
              return(cur_sd)
            })

            # which model is best ?
            idx_opt = which.min(lglik_vec)
            sd_opt  = err_vec[idx_opt]          
            fit_opt = models[[idx_opt]]
            
            # save coefficients
            pars$sd     = sd_opt
            psi0        = fit_opt$psi[2]
            if (is.null(psi0)) psi0 = NA
            names(psi0) = break.var
            pars$psi    = abs(psi0)
            
            if (class(fit_opt)[1] == 'segmented') {
              if (names(sd_opt) == 'right') {
                pars$coefficients <- fit_opt$coefficients[-length(fit_opt$coefficients)]
                str = paste('U1',break.var,sep='.')
                pars[[break.var]] = pars$coefficients[length(pars$coefficients)]
                names(pars$coefficients)[which(names(pars$coefficients) == str)] = break.var
                pars$coefficients[break.var] = 0
              }
              
              if (names(sd_opt) == 'left') {
                pars$coefficients <- fit_opt$coefficients[-length(fit_opt$coefficients)]
                str = paste('U1',break.var,sep='.')
                names(pars$coefficients)[which(names(pars$coefficients) == str)] = break.var
                pars[[break.var]] = 0
              }
              
              if (names(sd_opt) == 'dual') {
                pars$coefficients <- fit_opt$coefficients[-c(length(fit_opt$coefficients), length(fit_opt$coefficients)-1)]
                pars[[break.var]] = fit_opt$coefficients[length(fit_opt$coefficients)-1]
              }
            } else {
              
              if (names(sd_opt) == 'flat')  {
                pars$coefficients <- c(fit_opt$coefficients, 0)
                names(pars$coefficients) = c(names(fit_opt$coefficients), break.var)
                pars[[break.var]] = 0              
              }
              
              if (names(sd_opt) == 'slope')  {
                pars$coefficients <- fit_opt$coefficients
                pars[[break.var]] = 0              
              }              
            }
            
#             # perform tests for right slope
#             alpha        = 0.01
#             right.slope  = davies.test(fit.2, as.formula(paste('~',break.var)), k = 10)
#             right.slope  = right.slope$p.value <= alpha
#             
#             # perform test for left slope
#             tvals        = summary(fit_opt)$Ttable
#             left.slope   = tvals[break.var, 4] <= alpha
#             
#             # adjust coefficients
#             if (!right.slope) pars[[break.var]] = 0
#             if (!left.slope) pars$coefficients[break.var] = 0

            values        = unlist(pars)
            names(values) = c(names(pars$coefficients), "variance", names(pars$psi), break.var)
            object        = setpars(object, values)
            object@state_type   = names(sd_opt)
            object
          }
)

setMethod("setpars","breakpointFree",
          function(object, values, which="pars", prob=FALSE, ...) {
            npar <- npar(object)
            if(length(values)!=npar) stop("length of 'values' must be",npar)
            # determine whether parameters or fixed constraints are being set
            nms <- names(object@parameters$coefficients)
          
            if(length(values) == 0) return(object) # nothing to set; 
            switch(which,
                   "pars"= {
                     object@parameters$coefficients <- values[1:length(object@parameters$coefficients)] # matrix(values,ncol(object@x),byrow=TRUE) # this needs fixing!!!!                     
                     object@parameters$sd  <- values[(length(object@parameters$coefficients)+1)]
                     object@parameters$psi <- values[(length(object@parameters$coefficients)+2)]
                     object@parameters[[names(object@parameters$psi)]]<- values[(length(object@parameters$coefficients)+3)]                     
                   },
                   "fixed" = {
                     object@fixed <- as.logical(values)
                   }
            )
            names(object@parameters$coefficients) <- nms
            return(object)
          }
)

setMethod("getpars","breakpointFree",
          function(object,which="pars",...) {
            switch(which,
                   "pars" = {
                       parameters <- numeric()
                       parameters <- unlist(object@parameters)
                       pars <- parameters
                   },
                   "fixed" = {
                     pars <- object@fixed
                   }
            )
            return(pars)
          }
)

fitBreakpoint = function(data, fmla, K = 2, pstart = NULL, fixed = NULL, trstart = NULL, 
                         maxit = 20, tol = 5e-3) {
            
  # set up initializations
  if (is.null(trstart)) {
    trstart = matrix(rep(1, K^2), ncol = K) + diag(rep(K,K))
    trstart = trstart / rowSums(trstart)
  }

  breakpoint.var = all.vars(fmla)[2]
  params         = c('(Intercept)', all.vars(fmla)[-c(1)])
  break.val      = quantile(data[,breakpoint.var], probs = seq(0.25, 0.75, length.out = K))
  
  # set up response and transition for depmix fit
  rModels    = list()
  transition = list()
  for (k in 1:K) {
    
    if (is.null(pstart)) {
      pars = list(coefficients = rep(0.1, length(params)), sd = 1, psi = break.val[k], diff.slope = 0.1)
      pstart = unlist(pars)
      names(pstart) = NULL
      if (is.null(fixed)) fixed = as.logical(rep(0,length(pstart)))      
    }    
    
    bp = breakpointFree(formula = fmla, 
                        data    = data, 
                        pstart  = pstart, 
                        fixed   = fixed)
    rModels[[k]]    = list(bp)
    transition[[k]] = transInit(~1,nstates=K,data=data,pstart=trstart[k,])
  }            
  
  # set up initial distribution model
  instart = rep(1,K) / K
  inMod <- transInit(~1, ns=K, ps=instart, data=data.frame(1))
  
  mod <- makeDepmix(response=rModels, transition=transition, prior=inMod,stat=FALSE)
  
  fm  <- fit(mod,emc=em.control(rand=FALSE, maxit = maxit, tol = tol))
  
  return(fm)
}