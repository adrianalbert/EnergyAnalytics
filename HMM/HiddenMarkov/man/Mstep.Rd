\name{Mstep}
\alias{Mstep}
\alias{Mstep.beta}
\alias{Mstep.binom}
\alias{Mstep.exp}
\alias{Mstep.gamma}
\alias{Mstep.glm}
\alias{Mstep.lnorm}
\alias{Mstep.logis}
\alias{Mstep.norm}
\alias{Mstep.pois}

\title{M-Step of EM Algorithm for DTHMM}

\description{
Performs the \emph{maximisation} step of the EM algorithm for a \code{\link{dthmm}} process. This function is called by the \code{\link{BaumWelch}} function. The Baum-Welch algorithm used in the HMM literature is a version of the EM algorithm.
}
\usage{
Mstep.beta(x, cond, pm, pn, maxiter = 200)
Mstep.binom(x, cond, pm, pn)
Mstep.exp(x, cond, pm, pn)
Mstep.gamma(x, cond, pm, pn, maxiter = 200)
Mstep.glm(x, cond, pm, pn, family, link)
Mstep.lnorm(x, cond, pm, pn)
Mstep.logis(x, cond, pm, pn, maxiter = 200)
Mstep.norm(x, cond, pm, pn)
Mstep.pois(x, cond, pm, pn)
}

\arguments{
  \item{x}{is a vector of length \eqn{n} containing the observed process.}
  \item{cond}{is an object created by \code{\link{Estep}}.}
  \item{family}{character string, the GLM family, one of \code{"gaussian"}, \code{"poisson"}, \code{"Gamma"} or \code{"binomial"}.}
  \item{link}{character string, the link function. If \code{family == "Binomial"}, then one of \code{"logit"}, \code{"probit"} or \code{"cloglog"}; else one of \code{"identity"}, \code{"inverse"} or \code{"log"}.}
  \item{pm}{is a list object containing the current (Markov dependent) parameter estimates associated with the distribution of the observed process (see \code{\link{dthmm}}). These are only used as initial values if the algorithm within the \code{Mstep} is iterative.}
  \item{pn}{is a list object containing the observation dependent parameter values associated with the distribution of the observed process (see \code{\link{dthmm}}).}
  \item{maxiter}{maximum number of Newton-Raphson iterations.}
}

\details{
The functions \code{Mstep.beta}, \code{Mstep.binom}, \code{Mstep.exp}, \code{Mstep.gamma}, \code{Mstep.lnorm}, \code{Mstep.logis}, \code{Mstep.norm} and \code{Mstep.pois} perform the maximisation step for the Beta, Binomial, Exponential, Gamma, Log Normal, Logistic, Normal and Poisson distributions, respectively. Each function has the same argument list, even if specific arguments are redundant, because the functions are called from within other functions in a generic like manner. Specific notes for some follow.

\describe{
\item{\code{Mstep.beta}}{The \R functions for the \code{\link{Beta}} Distribution have arguments \code{shape1}, \code{shape2} and \code{ncp}. We only use \code{shape1} and \code{shape2}, i.e. \code{ncp} is assumed to be zero. Different combinations of \code{"shape1"} and \code{"shape2"} can be \dQuote{time} dependent (specified in \code{pn}) and Markov dependent (specified in \code{pm}). However, each should only be specified in one (see topic \code{\link{dthmm}}).
}

\item{\code{Mstep.binom}}{The \R functions for the \code{\link{Binomial}} Distribution have arguments \code{size} and \code{prob}. The \code{size} argument of the \code{\link{Binomial}} Distribution should always be specified in the \code{pn} argument (see topic \code{\link{dthmm}}).
}

\item{\code{Mstep.gamma}}{The \R functions for the \code{\link{GammaDist}} have arguments \code{shape}, \code{rate} and \code{scale}. Since \code{scale} is redundant, we only use \code{shape} and \code{rate}. Different combinations of \code{"shape"} and \code{"rate"} can be \dQuote{time} dependent (specified in \code{pn}) and Markov dependent (specified in \code{pm}). However, each should only be specified in one (see topic \code{\link{dthmm}}).
}

\item{\code{Mstep.lnorm}}{The \R functions for the \code{\link{Lognormal}} Distribution have arguments \code{meanlog} and \code{sdlog}. Different combinations of \code{"meanlog"} and \code{"sdlog"} can be \dQuote{time} dependent (specified in \code{pn}) and Markov dependent (specified in \code{pm}). However, each should only be specified in one (see topic \code{\link{dthmm}}).
}

\item{\code{Mstep.logis}}{The \R functions for the \code{\link{Logistic}} Distribution have arguments \code{location} and \code{scale}. Different combinations of \code{"location"} and \code{"scale"} can be \dQuote{time} dependent (specified in \code{pn}) and Markov dependent (specified in \code{pm}). However, each should only be specified in one (see topic \code{\link{dthmm}}).
}

\item{\code{Mstep.norm}}{The \R functions for the \code{\link{Normal}} Distribution have arguments \code{mean} and \code{sd}. Different combinations of \code{"mean"} and \code{"sd"} can be \dQuote{time} dependent (specified in \code{pn}) and Markov dependent (specified in \code{pm}). However, each should only be specified in one (see topic \code{\link{dthmm}}).
}

}
}

\section{Modifications and Extensions}{
The \pkg{HiddenMarkov} package calls the associated functions belonging to the specified probability distribution in a generic way. For example, if the argument \code{distn} in \code{\link{dthmm}} is \code{"xyz"}, it will expect to find functions \code{pxyz}, \code{dxyz}, and \code{Mstep.xyz}. And if simulations are to be performed, it will require \code{rxyz}. In this section we describe the required format for the distribution related functions \code{pxyz}, \code{dxyz}, and \code{rxyz}; and for the function \code{Mstep.xyz} required for the M-step in the EM algorithm.

Consider the examples below of distribution related functions and their arguments. Note that the probability functions all have a first argument of \code{q}, and the last two arguments are all the same, with the same default values. Similarly, the density functions have a first argument of \code{x}, and the last argument is the same, with the same defaults. The arguments in the middle are peculiar to the given distribution, \emph{one argument for each distribution parameter}. Note that the observed process \code{x} is \emph{univariate}.


\preformatted{
pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
pbeta(q, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
ppois(q, lambda, lower.tail = TRUE, log.p = FALSE)
pbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)

dnorm(x, mean = 0, sd = 1, log = FALSE)
dbeta(x, shape1, shape2, ncp = 0, log = FALSE)
dpois(x, lambda, log = FALSE)
dbinom(x, size, prob, log = FALSE)

rnorm(n, mean = 0, sd = 1)
rbeta(n, shape1, shape2, ncp = 0)
rpois(n, lambda)
rbinom(n, size, prob)

}

The functions \code{pxyz} (distribution function), \code{dxyz} (density) and \code{rxyz} (random number generator) must be consistent with the conventions used in the above examples. The software will deduce the distribution argument names from what is specified in \code{pm} and \code{pn}, and it will call these functions assuming that their argument list is consistent with those described above. The functions \code{pxyz} and \code{dxyz} are used in the forward and backward equations.
\cr

The functions \code{dxyz}, \code{pxyz} and \code{rxyz} \emph{must} also behave in the same vectorised way as \code{\link{dnorm}}. For example, if \code{x} is a vector, and \code{mean} and \code{sd} are scalars, then \code{\link{dnorm}(x, mean, sd)} calculates the density for each element in \code{x} using the scalar values of \code{mean} and \code{sd}; thus the returned value is the same length as \code{x}. Alternatively, if \code{x} is a scalar and \code{mean} and \code{sd} are vectors, both of the same length, then the returned value is the same length as \code{mean} and is the density of \code{x} evaluated at the corresponding pairs of values of \code{mean} and \code{sd}. The third possibility is that \code{x} and one of the distribution parameters, say \code{sd}, are vectors of the same length, and \code{mu} is a scalar. Then the returned vector will be of the same length as \code{x}, where the \eqn{i}th value is the density at \code{x[i]} with mean \code{mean} and standard deviation \code{sd[i]}. Note that the functions for the \code{\link{Multinomial}} distribution do not have this behaviour. Here the vector \code{x} contains the counts for \emph{one} multinomial experiment, so the vector is used to characterise the multivariate character of the random variable rather than multiple univariate realisations. Further, the distribution parameters (i.e. category probabilities) are characterised as one vector rather than a sequence of separate function arguments.
\cr

The other calculation, that is specific to the chosen distribution, is the maximisation in the M-step. If we have \code{distn="xyz"}, then there should be a function called \code{Mstep.xyz}. Further, it should have arguments \code{(x, cond, pm, pn)}; see for example \code{Mstep.norm}. The parameters that are estimated within this function are named in a consistent way with those that are defined within the \code{\link{dthmm}} arguments \code{pm} and \code{pn}. Notice that the estimates of \code{mean} and \code{sd} in \code{Mstep.norm} are weighted by \code{cond$u}. The calculations for \code{cond$u} are performed in the E-step, and utilise the distribution related functions \code{"dxyz"} and \code{"pxyz"}. The values of \code{cond$u} are essentially probabilities that the process belongs to the given Markov state, hence, when we calculate the distributional parameters (like \code{mu} and \code{sd} in \code{Mstep.norm}) we calculate weighted sums using the probabilities \code{cond$u}. This procedure can be shown to give the maximum likelihood estimates of \code{mu} and \code{sd}, and hence a similar weighting should be used for the distribution \code{"xyz"} (see Harte, 2006, for further mathematical detail). One needs to take a little more care when dealing with a distributions like the beta, where the cross derivatives of the log likelihood between the parameters, i.e. \eqn{\partial^2 \log L /(\partial \alpha_1 \partial \alpha_2)}{partial^2 log L /(partial alpha_1 partial alpha_2)} are non-zero. See \code{Mstep.beta} for further details.
\cr

Now consider a situation where we want to modify the way in which a normal distribution is fitted. Say we know the Markov dependent means, and we only want to estimate the standard deviations. Since both parameters are Markov dependent, they both need to be specified in the \code{pm} argument of \code{\link{dthmm}}. The estimation of the distribution specific parameters takes place in the M-step, in this case \code{\link{Mstep.norm}}. To achieve what we want, we need to modify this function. In this case it is relatively easy (see code in \dQuote{Examples} below). From the function \code{\link{Mstep.norm}}, take the code under the section \code{if (all(nms==c("mean", "sd")))}, i.e. both of the parameters are Markov dependent. However, replace the line where the mean is estimated to \code{mean <- pm$mean}, i.e. leave it as was initially specified. Unfortunately, one cannot easily modify the functions in a package namespace. The simple work-around here is to define a new distribution, say \code{"xyz"}, then define a new function with the above changes called \code{Mstep.xyz}. However, the distribution related functions are just the same as those for the normal distribution, hence, define them as follows:
\preformatted{
rxyz <- rnorm
dxyz <- dnorm
pxyz <- pnorm
qxyz <- qnorm
}
See the 2nd example below for full details.
}



\value{
A list object with the same structure as \code{pm} (see topic \code{\link{dthmm}}).
}

\seealso{
\code{\link{BaumWelch}}, \code{\link{Estep}}
}

\section{Note}{
The Mstep functions can be used to estimate the maximum likelihood parameters from a simple sample. See the example below where this is done for the logistic distribution.
}

\examples{
#    Fit logistic distribution to a simple single sample

#    Simulate data
n <- 20000
location <- -2
scale <- 1.5

x <- rlogis(n, location, scale)

#    give each datum equal weight
cond <- NULL
cond$u <- matrix(rep(1/n, n), ncol=1)

#    calculate maximum likelihood parameter estimates
#    start iterations at values used to simulate
print(Mstep.logis(x, cond,
                  pm=list(location=location,
                          scale=scale)))

#-----------------------------------------------------
#   Example with Gaussian Observations
#   Assume that both mean and sd are Markov dependent, but the means
#   are known and sd requires estimation (See "Modifications" above).
#   One way is to define a new distribution called "xyz", say.

Mstep.xyz <- function(x, cond, pm, pn){
    #   this function is a modified version of Mstep.norm
    #   here the mean is fixed to the values specified in pm$mean
    nms <- sort(names(pm))
    n <- length(x)
    m <- ncol(cond$u)
    if (all(nms==c("mean", "sd"))){
        mean <- pm$mean
        sd <- sqrt(apply((matrix(x, nrow=n, ncol=m) - 
                   matrix(mean,
                   nrow=n, ncol=m, byrow=TRUE))^2 * cond$u, MARGIN=2,
                   FUN=sum)/apply(cond$u, MARGIN=2, FUN=sum))
        return(list(mean=mean, sd=sd))
    }
}

#   define the distribution related functions for "xyz"
#   they are the same as those for the Normal distribution
rxyz <- rnorm
dxyz <- dnorm
pxyz <- pnorm
qxyz <- qnorm

Pi <- matrix(c(1/2, 1/2,   0,
               1/3, 1/3, 1/3,
                 0, 1/2, 1/2),
             byrow=TRUE, nrow=3)
p1 <- c(1, 6, 3)
p2 <- c(0.5, 1, 0.5)
n <- 1000

pm <- list(mean=p1, sd=p2)

x <- dthmm(NULL, Pi, c(0, 1, 0), "xyz", pm, discrete=FALSE)

x <- simulate(x, n, seed=5)

#    use above parameter values as initial values
y <- BaumWelch(x)

print(y$delta)
print(y$pm)
print(y$Pi)
}

\keyword{optimize}

