# below example recreates the same model as on the fit help page in a roundabout way
# there we had:
# mod1 <- depmix(list(rt~1,corr~1),data=speed,transition=~Pacc,nstates=2,
#   family=list(gaussian(),multinomial("identity")),ntimes=c(168,134,137))

data(speed)   

rModels <- list(
  list(
    GLMresponse(formula=rt~1,data=speed,family=gaussian()),
    GLMresponse(formula=corr~1,data=speed,family=multinomial("identity"))
  ),
  list(
    GLMresponse(formula=rt~1,data=speed,family=gaussian()),
    GLMresponse(formula=corr~1,data=speed,family=multinomial("identity"))
  )
)

transition <- list()
transition[[1]] <- transInit(~Pacc,nstates=2,data=speed)
transition[[2]] <- transInit(~Pacc,nstates=2,data=speed)

inMod <- transInit(~1,ns=2,data=data.frame(rep(1,3)),family=multinomial("identity"))
mod <- makeDepmix(response=rModels,transition=transition,prior=inMod,ntimes=c(168,134,137),stationary=FALSE)

set.seed(3)
fm1 <- fit(mod)
fm1
summary(fm1)


# generate data from two different multivariate normal distributions
m1 <- c(0,1)
sd1 <- matrix(c(1,0.7,.7,1),2,2)
m2 <- c(1,0)
sd2 <- matrix(c(2,.1,.1,1),2,2)
set.seed(2)
y1 <- mvrnorm(50,m1,sd1)
y2 <- mvrnorm(50,m2,sd2)
# this creates data with a single change point
y <- rbind(y1,y2)

# now use makeDepmix to create a depmix model for this bivariate normal timeseries
rModels <-  list()
rModels[[1]] <- list(MVNresponse(y~1))
rModels[[2]] <- list(MVNresponse(y~1))

trstart=c(0.9,0.1,0.1,0.9)

transition <- list()
transition[[1]] <- transInit(~1,nstates=2,data=data.frame(1),pstart=c(trstart[1:2]))
transition[[2]] <- transInit(~1,nstates=2,data=data.frame(1),pstart=c(trstart[3:4]))

instart=runif(2)
inMod <- transInit(~1,ns=2,ps=instart,data=data.frame(1))

mod <- makeDepmix(response=rModels,transition=transition,prior=inMod)

fm2 <- fit(mod,emc=em.control(random=FALSE))

# where is the switch point?
plot(as.ts(posterior(fm2)[,2]))


# in below example we add the exgaus distribution as a response model and fit
# this instead of the gaussian distribution to the rt slot of the speed data
# most of the actual computations for the exgaus distribution is done by calling
# functions from the gamlss family of packages; see their help pages for 
# interpretation of the mu, nu and sigma parameters that are fitted below

require(gamlss)
require(gamlss.dist)

data(speed)
rt <- speed$rt

# define a response class which only contains the standard slots, no additional slots
setClass("exgaus", contains="response")

# define a generic for the method defining the response class

setGeneric("exgaus", function(y, pstart = NULL, fixed = NULL, ...) standardGeneric("exgaus"))

# define the method that creates the response class

setMethod("exgaus", 
          signature(y="ANY"), 
          function(y,pstart=NULL,fixed=NULL, ...) {
            y <- matrix(y,length(y))
            x <- matrix(1)
            parameters <- list()
            npar <- 3
            if(is.null(fixed)) fixed <- as.logical(rep(0,npar))
            if(!is.null(pstart)) {
              if(length(pstart)!=npar) stop("length of 'pstart' must be ",npar)
              parameters$mu <- pstart[1]
              parameters$sigma <- log(pstart[2])
              parameters$nu <- log(pstart[3])
            }
            mod <- new("exgaus",parameters=parameters,fixed=fixed,x=x,y=y,npar=npar)
            mod
          }
)

setMethod("show","exgaus",
          function(object) {
            cat("Model of type exgaus (see ?gamlss for details) \n")
            cat("Parameters: \n")
            cat("mu: ", object@parameters$mu, "\n")
            cat("sigma: ", object@parameters$sigma, "\n")
            cat("nu: ", object@parameters$nu, "\n")
          }
)

setMethod("dens","exgaus",
          function(object,log=FALSE) {
            dexGAUS(object@y, mu = predict(object), sigma = exp(object@parameters$sigma), nu = exp(object@parameters$nu), log = log)
          }
)

setMethod("getpars","response",
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

setMethod("setpars","exgaus",
          function(object, values, which="pars", ...) {
            npar <- npar(object)
            if(length(values)!=npar) stop("length of 'values' must be",npar)
            # determine whether parameters or fixed constraints are being set
            nms <- names(object@parameters)
            switch(which,
                   "pars"= {
                     object@parameters$mu <- values[1]
                     object@parameters$sigma <- values[2]
                     object@parameters$nu <- values[3]
                   },
                   "fixed" = {
                     object@fixed <- as.logical(values)
                   }
            )
            names(object@parameters) <- nms
            return(object)
          }
)

setMethod("fit","exgaus",
          function(object,w) {
            if(missing(w)) w <- NULL
            y <- object@y
            fit <- gamlss(y~1,weights=w,family=exGAUS(),
                          control=gamlss.control(n.cyc=100,trace=FALSE),
                          mu.start=object@parameters$mu,
                          sigma.start=exp(object@parameters$sigma),
                          nu.start=exp(object@parameters$nu))
            pars <- c(fit$mu.coefficients,fit$sigma.coefficients,fit$nu.coefficients)
            object <- setpars(object,pars)
            object
          }
)

setMethod("predict","exgaus", 
          function(object) {
            ret <- object@parameters$mu
            return(ret)
          }
)

rModels <- list(
  list(
    exgaus(rt,pstart=c(5,.1,.1)),
    GLMresponse(formula=corr~1,data=speed,family=multinomial(),pstart=c(0.5,0.5))
  ),
  list(
    exgaus(rt,pstart=c(6,.1,.1)),
    GLMresponse(formula=corr~1,data=speed,family=multinomial(),pstart=c(.1,.9))
  )
)

trstart=c(0.9,0.1,0.1,0.9)
transition <- list()
transition[[1]] <- transInit(~Pacc,nstates=2,data=speed,pstart=c(trstart[1:2],0,0))
transition[[2]] <- transInit(~Pacc,nstates=2,data=speed,pstart=c(trstart[3:4],0,0))

instart=c(0.5,0.5)
inMod <- transInit(~1,ns=2,ps=instart,data=data.frame(rep(1,3)))

mod <- makeDepmix(response=rModels,transition=transition,prior=inMod,ntimes=c(168,134,137),stat=FALSE)

fm3 <- fit(mod,emc=em.control(rand=FALSE))
summary(fm3)