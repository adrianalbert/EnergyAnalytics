# responseresponseSelect.r
#
# Define breakpoint model for energy with temperature.
# 
# Adrian Albert
# Last modified: April 2013.
# -----------------------------------------------------------------------

library('segmented')

# define a response class which only contains the standard slots, no additional slots
setClass("responseSelect", 
         representation(formula    ="formula",
                        family     ="ANY",
                        state.type = 'character',
                        break.var  = 'character',
                        alpha      = 'numeric'),
         contains="response")

# define a generic for the method defining the response class
setGeneric("responseSelect",
           function(formula, state.type = NULL, break.var = NULL, alpha = NULL,
                    data = NULL, pstart = NULL, fixed = NULL, prob=TRUE, na.action="na.pass", ...) standardGeneric("responseSelect"))

# define the method that creates the response class
setMethod("responseSelect", 
          signature(formula = 'formula'), 
          function(formula, state.type = NULL, break.var = NULL, alpha = NULL,
                   data=NULL, pstart=NULL, fixed=NULL, prob=TRUE, na.action="na.pass", ...) {
            
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
                                    
            coefs.init = vector("numeric",length=ncol(x))
                        
            # set output parameters
            parameters <- list()
            constr <- NULL
            parameters$coefficients <- coefs.init
            parameters$sd <- 1
            constr <- list(
              parup = rep(Inf,length(parameters$coefficients)+1),
              parlow = c(rep(-Inf,length(parameters$coefficients)),0)
            )
            
            npar <- length(unlist(parameters))
            if(is.null(fixed)) fixed <- as.logical(rep(0,npar))
            if(!is.null(pstart)) {
              if(length(pstart)!=npar) stop("length of 'pstart' must be ",npar)
              parameters$coefficients[1:length(parameters$coefficients)] <- pstart[1:length(parameters$coefficients)]
              if(length(unlist(parameters))>length(parameters$coefficients)) {
                parameters$sd   <- as.numeric(pstart[(length(parameters$coefficients)+1)])
              }
            }                    
            
            if (is.null(state.type)) state.type = 'neutral'
            if (is.null(alpha)) alpha = '0.1'
            
            # mod <- new("exgaus",parameters=parameters,fixed=fixed,x=x,y=y,npar=npar)
            mod = new(Class = "responseSelect", formula = formula, 
                      state.type = state.type, break.var = break.var, alpha = alpha, 
                      parameters=parameters, fixed=fixed, x=x,y=y,npar=npar,constr=constr)
            # mod@state_type = ''
          }
)

setMethod("show","responseSelect",
          function(object) {
            cat("Model of type responseSelect formula: ", sep="")
            print(object@formula)            
            cat(sprintf("Type: %s; Dependent variable: %s\n",object@state.type, object@break.var))
            cat("Coefficients: \n")
            print(object@parameters$coefficients)
            cat("sd ",object@parameters$sd,"\n")
          }
)

setMethod("predict","responseSelect",
          function(object) {
            if (ncol(object@x) >= 2) {
              y.pred = object@x %*% object@parameters$coefficients
            } else 
              y.pred = object@x * object@parameters$coefficients
            return(y.pred)
          }
)

setMethod("dens","responseSelect",
          function(object,log=FALSE) {
            dnorm(x=object@y,mean=predict(object),sd=object@parameters$sd,log=log)
          }
)

setGeneric("logDens",function(object,...) standardGeneric("logDens"))
setMethod("logDens","responseSelect",
          function(object) {
            dnorm(x=object@y,mean=predict(object),sd=object@parameters$sd,log=TRUE)
          }
)

setMethod("logLik","responseSelect",
          function(object) {
            sum(logDens(object))
          }
)

setMethod("simulate",signature(object="responseSelect"),
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

setMethod("fit","responseSelect",
          function(object, w, maxit = 50, tol = 1e-2) {
            
            if(missing(w) | length(which(w>0))==0) w <- NULL
            nas <- is.na(rowSums(object@y))
            pars <- object@parameters
            break.var = object@break.var
            
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

            # compute best model
            models = list(flat = fit.1, slope = fit.2) 
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
            
            if (names(sd_opt) == 'flat')  {
              pars$coefficients <- c(fit_opt$coefficients, 0)
              names(pars$coefficients) = c(names(fit_opt$coefficients), break.var)
              object@state.type = 'neutral'
            }
            
            if (names(sd_opt) == 'slope')  {
              pars$coefficients <- fit_opt$coefficients
              # perform tests for slope
              alpha        = object@alpha
              pvals        = summary(fit_opt)$coefficients[,4]
              has.slope    = pvals[break.var] <= alpha
              
              # adjust coefficients
              if (!has.slope) {
                pars$coefficients[break.var] = 0
                object@state.type = 'neutral'
              } else {
                if (pars$coefficients[break.var] < 0) object@state.type = 'heating' else object@state.type = 'cooling'
              }              
            }              
            
            values        = unlist(pars)
            names(values) = c(names(pars$coefficients), "variance")
            object        = setpars(object, values)
            object
          }
)

setMethod("setpars","responseSelect",
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
                   },
                   "fixed" = {
                     object@fixed <- as.logical(values)
                   }
            )
            names(object@parameters$coefficients) <- nms
            return(object)
          }
)

setMethod("getpars","responseSelect",
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

makeModelSelect = function(data, fmla.resp, fmla.tran, 
                         K = 2, alpha = 0.001) {
            
  breakpoint.var = all.vars(fmla.resp)[2]
  params.resp    = all.vars(fmla.resp)[-1]
  params.resp    = c('(Intercept)', unlist(sapply(params.resp, function(p) paste(p,levels(data[,p])[-1],sep=''))))
  params.tran    = all.vars(fmla.tran)
  params.tran    = c('(Intercept)', unlist(sapply(params.tran, function(p) paste(p,levels(data[,p])[-1],sep=''))))
  
  # set up response and transition for depmix fit
  rModels    = list()
  transition = list()
  for (k in 1:K) {
    
    # initialize response model
    pars = list(coefficients = runif(length(params.resp)), sd = 1)
    pstart = unlist(pars)
    names(pstart) = NULL
    fixed = as.logical(rep(0,length(pstart)))      
    
    bp = responseSelect(formula   = fmla.resp, 
                        data      = data, 
                        pstart    = pstart, 
                        fixed     = fixed,
                        break.var = breakpoint.var,
                        alpha     = alpha,
                        state.type= 'neutral')
    rModels[[k]]    = list(bp)
    
    
    # initialize transition model
    nvars.tran = length(params.tran)+1
    trstart = NULL# matrix(runif(K*nvars.tran), ncol = K)
    transition[[k]] = transInit(fmla.tran, nstates=K, data=data, 
                                pstart=trstart, fixed = as.logical(c(1,rep(0,(K-1)*nvars.tran))))
  }            
  
  # set up initial distribution model
  instart = rep(1,K) / K
  inMod <- transInit(~1, ns=K, ps=instart)
  
  # define and fit model
  mod <- makeDepmix(response=rModels, transition=transition, prior=inMod, stat=FALSE)  
  
  return(mod)
}