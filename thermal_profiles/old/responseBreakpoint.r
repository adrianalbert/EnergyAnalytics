# responseBreakpoint.r
#
# Define breakpoint model for energy with temperature.
# 
# Adrian Albert
# Last modified: April 2013.
# -----------------------------------------------------------------------

library('segmented')

# define a response class which only contains the standard slots, no additional slots
setClass("breakpoint", 
         representation(formula="formula",
                        family="ANY"),
         contains="response")

# define a generic for the method defining the response class
setGeneric("breakpoint",
           function(formula, breakpoint.var = NULL, data = NULL, pstart = NULL, fixed = NULL, prob=TRUE, na.action="na.pass", ...) standardGeneric("breakpoint"))

# define the method that creates the response class
setMethod("breakpoint", 
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
            fixed[length(fixed)] = TRUE
            if(!is.null(pstart)) {
              if(length(pstart)!=npar) stop("length of 'pstart' must be ",npar)
              parameters$coefficients[1:length(parameters$coefficients)] <- pstart[1:length(parameters$coefficients)]
              if(length(unlist(parameters))>length(parameters$coefficients)) {
                parameters$sd   <- as.numeric(pstart[(length(parameters$coefficients)+1)])
                parameters$psi  <- as.numeric(pstart[(length(parameters$coefficients)+2)])
                parameters[[breakpoint.var]] <- as.numeric(pstart[(length(parameters$coefficients)+3)])
              }
            }                    
            
            # mod <- new("exgaus",parameters=parameters,fixed=fixed,x=x,y=y,npar=npar)
            mod = new(Class = "breakpoint", formula = formula, parameters=parameters, fixed=fixed, x=x,y=y,npar=npar,constr=constr)
            mod
          }
)

setMethod("show","breakpoint",
          function(object) {
            cat("Model of type breakpoint (see ?segmented for details) formula: ", sep="")
            print(object@formula)
            cat("Coefficients: \n")
            print(object@parameters$coefficients)
            cat("sd ",object@parameters$sd,"\n")
            cat('Break point\n')
            print(object@parameters$psi)
            cat('Difference in slope\n')
            print(object@parameters[[4]])
          }
)

setMethod("predict","breakpoint",
          function(object) {
            break.var = names(object@parameters)[4]
            idx    = which(colnames(object@x) == break.var)
            if (ncol(object@x) > 2) {
              y.pred = object@x %*% object@parameters$coefficients
            } else 
              y.pred = object@x %*% object@parameters$coefficients
            indict = (object@x[,idx] - object@parameters$psi) > 0
            y.adj  = (object@x[,idx] - object@parameters$psi) * object@parameters[[break.var]] * indict
            y.pred = y.pred + y.adj
            return(y.pred)
          }
)

setMethod("dens","breakpoint",
          function(object,log=FALSE) {
            dnorm(x=object@y,mean=predict(object),sd=object@parameters$sd,log=log)
          }
)

setGeneric("logDens",function(object,...) standardGeneric("logDens"))
setMethod("logDens","breakpoint",
          function(object) {
            dnorm(x=object@y,mean=predict(object),sd=object@parameters$sd,log=TRUE)
          }
)

setMethod("logLik","breakpoint",
          function(object) {
            sum(logDens(object))
          }
)

setMethod("simulate",signature(object="breakpoint"),
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
setMethod("fit","breakpoint",
          function(object, w) {
            if(missing(w)) w <- NULL
            nas <- is.na(rowSums(object@y))
            pars <- object@parameters
            break.var = names(pars)[4]
            
            # first prepare data for regular lm model
            data          = as.data.frame(cbind(object@x[!nas,], object@y[!nas,]))
            names(data)   = c(colnames(object@x), colnames(object@y))            
            idx_intercept = which(names(data) == '(Intercept)')
            if (length(idx_intercept)>0) data = data[,-idx_intercept]

            # is it one-sided or two-sided breakpoint model?
#             predictors    = intersect(names(data), setdiff(colnames(object@x), break.var))
            predictors    = intersect(names(data), colnames(object@x))
            
            # fit regular lm model
            if (length(predictors)>0) {
              fmla_lm     = as.formula(paste(colnames(object@y), '~', paste(predictors, collapse = '+')))
            } else 
              fmla_lm     = as.formula(paste(colnames(object@y), '~ 1'))

            if(!is.null(w)) {              
              fit <- lm(fmla_lm, data = data, weights = w[!nas])
            } else {
              fit <- lm(fmla_lm, data = data)
            }
            
            # then compute segmented model
            fmla = as.formula(paste('~',break.var))
#             fit.seg<-segmented(fit, seg.Z = fmla, psi = pars$psi)            
            fit.seg <- try(segmented(fit, seg.Z = fmla, psi = pars$psi))
            
            # save coefficients
            if (class(fit.seg) == 'try-error') {
              pars$coefficients <- fit$coefficients
              if(!is.null(w)) {
                pars$sd <- sqrt(sum(w[!nas]*fit$residuals^2/sum(w[!nas])))
              } else {
                pars$sd <- sd(fit$residuals)
              }
              # pars$psi  = 0
              pars[[4]] = 0              
            } else {
              pars$coefficients <- fit.seg$coefficients[-c(length(fit.seg$coefficients)-1, length(fit.seg$coefficients))]
              if(!is.null(w)) {
                pars$sd <- sqrt(sum(w[!nas]*fit.seg$residuals^2/sum(w[!nas])))
              } else {
                pars$sd <- sd(fit.seg$residuals)
              }
              pars$psi  = fit.seg$psi[2]
              pars[[4]] = fit.seg$coefficients[length(fit.seg$coefficients)-1]
            }
            
            object <- setpars(object,unlist(pars))
            object
          }
)

setMethod("setpars","breakpoint",
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
                     object@parameters[[4]]<- values[(length(object@parameters$coefficients)+3)]
                   },
                   "fixed" = {
                     object@fixed <- as.logical(values)
                   }
            )
            names(object@parameters$coefficients) <- nms
            return(object)
          }
)

setMethod("getpars","breakpoint",
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
