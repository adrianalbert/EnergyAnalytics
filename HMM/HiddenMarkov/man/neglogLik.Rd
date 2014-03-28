\name{neglogLik}
\alias{neglogLik}

\title{Negative Log-Likelihood}
\description{
Calculates the log-likelihood multiplied by negative one. It is in a format that can be used with the functions \code{\link[stats]{nlm}} and \code{\link[stats]{optim}}, providing an alternative to the \code{\link{BaumWelch}} algorithm for maximum likelihood parameter estimation.
}

\usage{
neglogLik(params, object, pmap)
}

\arguments{
  \item{params}{a vector of revised parameter values.}
  \item{object}{an object of class \code{"\link{dthmm}"}, \code{"\link{mmglm0}"}, or \code{"\link{mmpp}"}.}
  \item{pmap}{a user provided function mapping the revised (or restricted) parameter values \code{p} into the appropriate locations in \code{object}.}
}

\details{
This function is in a format that can be used with the two functions \code{\link[stats]{nlm}} and \code{\link[stats]{optim}} (see Examples below). This provides alternative methods of estimating the maximum likelihood parameter values, to that of the EM algorithm provided by \code{\link{BaumWelch}}, including Newton type methods and grid searches. It can also be used to restrict estimation to a subset of parameters.

The EM algorithm is relatively stable when starting from poor initial values but convergence is very slow in close proximity to the solution. Newton type methods are very sensitive to initial conditions but converge much more quickly in close proximity to the solution. This suggests initially using the EM algorithm and then switching to Newton type methods (see Examples below).

The maximisation of the model likelihood function can be restricted to be over a subset of the model parameters. Other parameters will then be fixed at the values stored in the model \code{object}. Let \eqn{\Theta}{Theta} denote the model parameter space, and let \eqn{\Psi}{Psi} denote the parameter sub-space (\eqn{\Psi \subseteq \Theta}{Psi subseteq Theta}) over which the likelihood function is to be maximised. The argument \code{params} contains values in \eqn{\Psi}{Psi}, and \code{pmap} is assigned a function that maps these values into the model parameter space \eqn{\Theta}{Theta}. See \dQuote{Examples} below.

The mapping function assigned to \code{pmap} can also be made to impose restrictions on the domain of the parameter space \eqn{\Psi}{Psi} so that the minimiser cannot jump to values such that \eqn{\Psi \not\subseteq \Theta}{Psi notsubseteq Theta}. For example, if a particular parameter must be positive, one can work with a transformed parameter that can take any value on the real line, with the model parameter being the exponential of this transformed parameter. Similarly a modified logit like transform can be used to ensure that parameter values remain within a fixed interval with finite boundaries. Examples of these situations can be found in the  \dQuote{Examples} below.

Some functions are provided in the topic \code{\link{Transform-Parameters}} that may provide useful components within the user provided function assigned to \code{pmap}.
}

\value{
Value of the log-likelihood times negative one.
}

\seealso{\code{\link[stats]{nlm}}, \code{\link[stats]{optim}}, \code{\link{Transform-Parameters}}, \code{\link{BaumWelch}}}

\examples{
#   Example where the Markov chain is assumed to be stationary

Pi <- matrix(c(0.8, 0.1, 0.1,
               0.1, 0.6, 0.3,
               0.2, 0.3, 0.5),
             byrow=TRUE, nrow=3)

#   start simulation in state 2
delta <- c(0, 1, 0)

x <- dthmm(NULL, Pi, delta, "exp", list(rate=c(5, 2, 0.2)), nonstat=FALSE)
x <- simulate(x, nsim=5000, seed=5)

#   Approximate solution using BaumWelch
x1 <- BaumWelch(x, control=bwcontrol(maxiter=10, tol=1e-5))


#   Exact solution using nlm

allmap <- function(y, p){
    #    maps vector back to delta, Pi and rate
    m <- sqrt(length(p))
    y$Pi <- vector2Pi(p[1:(m*(m-1))])
    y$pm$rate <- exp(p[(m^2-m+1):(m^2)])
    y$delta <- compdelta(Pi)
    return(y)
}

p <- c(Pi2vector(x$Pi), log(x$pm$rate))
#   Increase iterlim below to achieve convergence
#   Has been restricted to minimise time of package checks
z <- nlm(neglogLik, p, object=x, pmap=allmap,
         print.level=2, gradtol=0.000001, iterlim=2)
x2 <- allmap(x, z$estimate)


#   compare parameter estimates
print(summary(x))
print(summary(x1))
print(summary(x2))


#--------------------------------------------------------
#   Estimate only the off diagonal elements in the matrix Pi
#   Hold all others as in the simulation

#   This function maps the changeable parameters into the
#   dthmm object - done within the function neglogLik
#   The logit-like transform removes boundaries

Pi <- matrix(c(0.8, 0.1, 0.1,
               0.1, 0.6, 0.3,
               0.2, 0.3, 0.5),
             byrow=TRUE, nrow=3)

delta <- c(0, 1, 0)

x <- dthmm(NULL, Pi, delta, "exp", list(rate=c(5, 3, 1)))
x <- simulate(x, nsim=5000, seed=5)

offdiagmap <- function(y, p){
    #   rows must sum to one
    invlogit <- function(eta)
        exp(eta)/(1+exp(eta))
    y$Pi[1,2] <- (1-y$Pi[1,1])*invlogit(p[1])
    y$Pi[1,3] <- 1-y$Pi[1,1]-y$Pi[1,2]
    y$Pi[2,1] <- (1-y$Pi[2,2])*invlogit(p[2])
    y$Pi[2,3] <- 1-y$Pi[2,1]-y$Pi[2,2]
    y$Pi[3,1] <- (1-y$Pi[3,3])*invlogit(p[3])
    y$Pi[3,2] <- 1-y$Pi[3,1]-y$Pi[3,3]
    return(y)
}

z <- nlm(neglogLik, c(0, 0, 0), object=x, pmap=offdiagmap,
         print.level=2, gradtol=0.000001)

#    x1 contains revised parameter estimates
x1 <- offdiagmap(x, z$estimate)

#    print revised values of Pi
print(x1$Pi)

#    print log-likelihood using original and revised values
print(logLik(x))
print(logLik(x1))

#--------------------------------------------------------
#   Fully estimate both Q and lambda for an MMPP Process

Q <- matrix(c(-8,  5,  3,
               1, -4,  3,
               2,  5, -7),
            byrow=TRUE, nrow=3)/25
lambda <- c(5, 3, 1)
delta <- c(0, 1, 0)

#    simulate some data
x <- mmpp(NULL, Q, delta, lambda)
x <- simulate(x, nsim=5000, seed=5)

allmap <- function(y, p){
    #    maps vector back to Pi and rate
    m <- sqrt(length(p))
    y$Q <- vector2Q(p[1:(m*(m-1))])
    y$lambda <- exp(p[(m^2-m+1):(m^2)])
    return(y)
}

#    Start by using the EM algorithm
x1 <- BaumWelch(x, control=bwcontrol(maxiter=10, tol=0.01))

#    use above as initial values for the nlm function
#    map parameters to a single vector, fixed delta
p <- c(Q2vector(x1$Q), log(x1$lambda))

#    Complete estimation using nlm
#    Increase iterlim below to achieve convergence
#    Has been restricted to minimise time of package checks
z <- nlm(neglogLik, p, object=x, pmap=allmap,
         print.level=2, gradtol=0.000001, iterlim=5)

#    mmpp object with estimated parameter values from nlm
x2 <- allmap(x, z$estimate)

#    compare log-likelihoods
print(logLik(x))
print(logLik(x1))
print(logLik(x2))

#   print final parameter estimates
print(summary(x2))
}

\keyword{optimize}
