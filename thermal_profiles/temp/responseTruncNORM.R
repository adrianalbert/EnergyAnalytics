# responseTruncNORM.r
#
# Define response as truncated normal.
# 
# Adrian Albert
# Last modified: October 2013.
# -----------------------------------------------------------------------

library('truncnorm') # truncated normal distribution
library('VGAM')      # censored regression

setClass("TruncNORMresponse",contains="GLMresponse")

# method 'fit'
# use: in EM (M step)
# returns: (fitted) response with (new) estimates of parameters
# The M step here is a Tobit (censored) regression, since the dependent variable is always >0.

setMethod("fit","TruncNORMresponse",
	function(object,w) {
		if(missing(w)) w <- NULL
    nas <- is.na(rowSums(object@y))
	  pars <- object@parameters
		dat <- cbind(object@x[!nas,], object@y[!nas,])
		fml <- as.formula(paste(names(object@y[!nas,]), '~', paste(names(object@y[!nas,]), collapse='+')))		
		if(!is.null(w)) {
		  fit <- vglm(fml, tobit(Lower = 0), data = dat, weights = w[!nas])
			# fit <- lm.wfit(x=as.matrix(object@x[!nas,]),y=as.matrix(object@y[!nas,]),w=w[!nas])
		} else {
# 			fit <- lm.fit(x=as.matrix(object@x[!nas,]),y=as.matrix(object@y[!nas,]))
		  fit <- vglm(fml, tobit(Lower = 0), data = dat)
		}
		pars$coefficients <- fit$coefficients
		if(!is.null(w)) {
			pars$sd <- sqrt(sum(w[!nas]*fit$residuals^2/sum(w[!nas])))
		} else {
			pars$sd <- sd(fit$residuals)
		}
		object <- setpars(object,unlist(pars))
		object
	}
)

setMethod("logDens","TruncNORMresponse",
	function(object) {
	  dtruncnorm(object@y, a = 0, mean=predict(object),sd=object@parameters$sd,log=TRUE)
	}
)

setMethod("dens","TruncNORMresponse",
	function(object,log=FALSE) {
		dtruncnorm(x=object@y, a = 0, mean=predict(object), sd=object@parameters$sd, log=log)
	}
)

setMethod("predict","TruncNORMresponse",
	function(object) {
		object@x%*%object@parameters$coefficients
	}
)

setMethod("simulate",signature(object="TruncNORMresponse"),
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
    response <- rtruncnorm(nt*nsim, a = 0, mean=mu,sd=sd)
    #if(nsim > 1) response <- matrix(response,ncol=nsim)
    response <- as.matrix(response)
    return(response)
  }
)

# Create depmixS4 model to use for fitting custom distributions
makeModelTruncNORM = function(data, fmla.resp, fmla.tran, K = 2) {
  
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
    
    bp = new("TruncNORMresponse", data = data, formula = fmla.resp)
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