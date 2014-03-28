\name{mmglm}
\alias{mmglm0}
\alias{mmglm1}
\alias{mmglmlong1}

\title{Markov Modulated GLM Object}

\description{
These functions create Markov modulated generalised linear model objects. \emph{\bold{These functions are in development and may change}}, see \dQuote{Under Development} below.
}

\usage{
mmglm0(x, Pi, delta, family, link, beta, glmformula = formula(y~x1),
       sigma = NA, nonstat = TRUE, msg = TRUE)
mmglm1(y, Pi, delta, glmfamily, beta, Xdesign,
       sigma = NA, nonstat = TRUE, size = NA, msg = TRUE)
mmglmlong1(y, Pi, delta, glmfamily, beta, Xdesign, longitude, 
           sigma = NA, nonstat = TRUE, size = NA, msg = TRUE)
}

\arguments{
  \item{x}{a dataframe containing the observed variable (i.e. the response variable in the generalised linear model) and the covariate. The function \code{mmglm0} requires that the response variable be named \code{y} and the covariate \code{x1}.  Alternatively, \code{x} could be specified as \code{NULL}, meaning that the data will be added later (e.g. simulated). See Details below for the binomial case. The functions \code{mmglm1} and \code{mmglmlong1} do not have these naming restrictions.}
  \item{y}{numeric vector, response variable. In the case of binomial, it is the number of successes (see argument \code{size}).}
  \item{Pi}{is the \eqn{m \times m}{m*m} transition probability matrix of the hidden Markov chain.}
  \item{delta}{is the marginal probability distribution of the \eqn{m} hidden states at the first time point.}
  \item{family}{character string, the GLM family, one of \code{"gaussian"}, \code{"poisson"}, \code{"Gamma"} or \code{"binomial"}.}
  \item{link}{character string, the link function. If \code{family == "binomial"}, then one of \code{"logit"}, \code{"probit"} or \code{"cloglog"}; else one of \code{"identity"}, \code{"inverse"} or \code{"log"}.}
  \item{glmfamily}{a \code{\link[stats]{family}} object defining the glm family and link function. It is currently restricted to Gaussian, Poisson, Binomial or Gamma models with the standard link functions provided by \code{\link[stats]{glm}}.}
  \item{Xdesign}{a \eqn{nN \times p}{(nN)*p} design matrix, where \eqn{p} is the number of parameters in the linear predictor, \eqn{N} is the number of subjects (\eqn{N=1} in \code{mmglm1}), and \eqn{n} is the number of observations for each subject (\emph{assumed to be the same}).}
  \item{beta}{a \eqn{p \times m}{p*m} matrix containing parameter values, used as initial values during estimation. In the case of the simple regression model of \code{mmglm0}, \eqn{p=2}. In the case of \code{mmglm1} and \code{mmglmlong1}, \eqn{p} is the number of columns of \code{Xdesign}.}
  \item{glmformula}{the only model formula for \code{mmglm0} is \code{y~x1}. Note that the functions \code{mmglm1} and \code{mmglmlong1} do not have this restriction, however, in those cases, the model formula is currently implicitly defined through \code{Xdesign}.}
  \item{sigma}{if \code{family == "gaussian"}, then it is the variance; if \code{family == "Gamma"}, then it is \code{1/sqrt(shape)}. It is of length \eqn{m} for each Markov state.}
  \item{nonstat}{is logical, \code{TRUE} if the homogeneous Markov chain is assumed to be non-stationary, default.}
  \item{longitude}{a vector the same length as \code{y} identifying the subject for each observation. The observations must be grouped by subject, and ordered by \dQuote{time} within subject.}
  \item{size}{is number of Bernoulli trials in each observation when the glm \code{\link[stats]{family}} is binomial. It is the same length as \code{y}.}
  \item{msg}{is logical, suppress messages about developmental status.}
}

\value{
A \code{\link[base]{list}} object with class \code{"mmglm0"}, containing the above arguments as named components.
}

\details{
This family of models is similar in nature to those of the class \code{\link{dthmm}}, in that both classes have the distribution of the observed variable being \dQuote{modulated} by the changing hidden Markov state. They differ slightly in the mechanism. This family assumes that the mean of the observation distribution can be expressed as a linear model of other known variables, but it is the parameters in the linear predictor that are being modulated by the hidden Markov process, thus causing the changes in the observed means. The linear model is assumed to be a generalised linear model as described by McCullagh & Nelder (1989).

The function \code{mmglm0} is a very simple trivial case where the linear predictor is of the form \eqn{\beta_0 + \beta_1 x_1}{beta0 + beta1*x1}. The version \code{mmglm1} does not have this limitation. The model formula for \code{mmglm1} is defined implicitly through the structure of the specified design matrix. The model \code{mmglmlong1} is similar to \code{mmglm1} but can be applied to longitudinal observations. Models of the form given by \code{mmglm1} are assumed to have one time series, and from a theoretical perspective, one would be interested in the asymptotic properties of the parameter estimates as the series length gets very large. In the longitudinal case (\code{mmglmlong1}), the series of observations per individual is probably very small (\eqn{<10}), and hence interest is in the asymptotic properties as the number of individuals becomes large. Note that in the longitudinal case, the number of observations per individual is assumed to be the same. The responses are assumed to be conditionally independent given the value of the Markov chain and the explanatory variables in the linear predictor.

If \code{family == "binomial"} then the response variable \code{y} is interpreted as the number of successes. The dataframe \code{x} must also contain a variable called \code{size} being the number of Bernoulli trials. This is different to the format used by the function \code{\link{glm}} where \code{y} would be a matrix with two columns containing the number of successes and failures, respectively. The different format here allows one to specify the number of Bernoulli trials \emph{only} so that the number of successes or failures can be simulated later.

When the density function of the response variable is from the exponential family (Charnes et al, 1976, Eq. 2.1), the likelihood function (Charnes et al, 1976, Eq. 2.4) can be maximised by using iterative weighted least squares (Charnes et al, 1976, Eq. 1.1 and 1.2). This is the method used by the \R function \code{\link{glm}}. In this Markov modulated version of the model, the third term of the complete data log-likelihood, as given in Harte (2006, Sec. 2.3), needs to be maximised. This is simply the sum of the individual log-likelihood contributions of the response variable weighted by the Markov state probabilities calculated in the E-step. This can also be maximised using iterative least squares by passing these additional weights (Markov state probabilities) into the \code{\link{glm}} function.
}

\section{Under Development}{
These functions are still being developed. In previous releases of the package (\eqn{< 1.3}), there was only one function called \code{mmglm}. This has been renamed to \code{mmglm0}. The most recent version is \code{mmglm1} along with \code{mmglmlong1} which has flexibility to include longitudinal data. Further development versions will be numbered sequentially. The name \code{mmglm} has been reserved for the final stable version, at which point the numbered versions will become deprecated.
}

\references{
Cited references are listed on the \link{HiddenMarkov} manual page.
}

\examples{
#--------------------------------------------------------
#     Gaussian with identity link function
#         using mmglm0

delta <- c(0,1)

Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

beta <- matrix(c(0.1, -0.1,
                 1.0,  5.0),
               byrow=TRUE, nrow=2)

x <- mmglm0(NULL, Pi, delta, family="gaussian", link="identity",
            beta=beta, sigma=c(1, 2))

n <- 1000
x <- simulate(x, nsim=n, seed=10)

#   Increase maxiter below to achieve convergence
#   Has been restricted to minimise time of package checks
y <- BaumWelch(x, bwcontrol(maxiter=2))

w <- hist(residuals(y))
z <- seq(-3, 3, 0.01)
points(z, dnorm(z)*n*(w$breaks[2]-w$breaks[1]), col="red", type="l")
box()

print(summary(y))
print(logLik(y))


#--------------------------------------------------------
#    Gaussian with log link function
#         using mmglm1

n <- 1000

#   the range of x needs changing according to the glmfamily
x <- seq(-0.9, 1.5, length.out=n)

colour <- c("blue", "green", "red")
colnum <- rep(1:3, n/3+1)[1:n] - 1

data <- data.frame(x=x, colour=colour[colnum+1])

#   will simulate response variable, not required in formula
#   design matrix only depends on RHS of formula
glmformula <- formula( ~ x + I(x^2) + colour)
glmfamily <- gaussian(link="log")
Xdesign <- model.matrix(glmformula, data=data)

# --- Parameter Values and Simulation ---

Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

delta <- c(1, 0)

sd <- c(1.2, 1)

beta <- matrix(c(-1, -1.2,
                 -2, -1.8,
                  3,  2.8,
                  1,  0.8, 
                  2,  2.2), 
               ncol=ncol(Pi), nrow=ncol(Xdesign), byrow=TRUE)

y <- mmglm1(NULL, Pi, delta, glmfamily, beta, Xdesign, sigma=sd)

y <- simulate(y, seed=5)

# --- Estimation ---

#   Increase maxiter below to achieve convergence
#   Has been restricted to minimise time of package checks
tmp <- BaumWelch(y, bwcontrol(posdiff=FALSE, maxiter=2))
print(summary(tmp))


#-------------------------------------------------
#    Binomial with logit link function
#         using mmglm1

#   n = series length
n <- 1000

#   the range of x need changing according to the glmfamily
x <- seq(-1, 1.5, length.out=n)

colour <- c("blue", "green", "red")
colnum <- rep(1:3, n/3+1)[1:n] - 1

data <- data.frame(x=x, colour=colour[colnum+1])

glmformula <- formula( ~ x + I(x^2) + colour)
glmfamily <- binomial(link="logit")
Xdesign <- model.matrix(glmformula, data=data)

# --- Parameter Values and Simulation ---

Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

delta <- c(1, 0)

beta <- matrix(c(-1, -1.2,
                 -2, -1.8,
                  3,  2.8,
                  1,  0.8, 
                  2,  2.2), 
               ncol=ncol(Pi), nrow=ncol(Xdesign), byrow=TRUE)

y <- mmglm1(NULL, Pi, delta, glmfamily, beta, Xdesign, sigma=sd,
            size=rep(100, n))

#   each element of y$y is the number of successes in 100 Bernoulli trials
y <- simulate(y, seed=5)


# --- Estimation ---

#   Increase maxiter below to achieve convergence
#   Has been restricted to minimise time of package checks
tmp <- BaumWelch(y, bwcontrol(posdiff=FALSE, maxiter=2))
print(summary(tmp))


#-------------------------------------------------
#    Gaussian with log link function, longitudinal data
#         using mmglmlong1

#   n = series length for each subject
#   N = number of subjects
n <- 5
N <- 1000

#   the range of x need changing according to the glmfamily
x <- seq(-0.9, 1.5, length.out=n)

colour <- c("blue", "green", "red")
colnum <- rep(1:3, n/3+1)[1:n] - 1

data <- data.frame(x=x, colour=colour[colnum+1])

#   will simulate response variable, not required in formula
#   design matrix only depends on RHS of formula
glmformula <- formula( ~ x + I(x^2) + colour)
glmfamily <- gaussian(link="log")
Xdesign0 <- model.matrix(glmformula, data=data)

#    multiple subjects
Xdesign <- NULL
for (i in 1:N) Xdesign <- rbind(Xdesign, Xdesign0)

# --- Parameter Values and Simulation ---

Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

delta <- c(0.5, 0.5)

sd <- c(1.2, 1)

beta <- matrix(c(-1, -1.2,
                 -2, -1.8,
                  3,  2.8,
                  1,  0.8, 
                  2,  2.2), 
               ncol=ncol(Pi), nrow=ncol(Xdesign), byrow=TRUE)

y <- mmglmlong1(NULL, Pi, delta, glmfamily, beta, Xdesign, sigma=sd,
                longitude=rep(1:N, each=n))

y <- simulate(y, seed=5)

# --- Estimation ---

#    Note: the "Not run" blocks below are not run during package checks
#    as the makeSOCKcluster definition is specific to my network,
#    modify accordingly if you want parallel processing.

cl <- NULL
\dontrun{
if (require(snow)){
    cl <- makeSOCKcluster(c("localhost", "horoeka.localdomain", 
                            "horoeka.localdomain", "localhost"))
}}

#   Increase maxiter below to achieve convergence
#   Has been restricted to minimise time of package checks
tmp <- BaumWelch(y, bwcontrol(posdiff=FALSE, maxiter=2),
                 SNOWcluster=cl)

\dontrun{
if (!is.null(cl)){
    stopCluster(cl)
    rm(cl)
}}

print(summary(tmp))


#-------------------------------------------------
#    Binomial with logit link function, longitudinal data
#         using mmglmlong1

#   n = series length for each subject
#   N = number of subjects
n <- 10
N <- 100

#   the range of x need changing according to the glmfamily
x <- seq(-1, 1.5, length.out=n)

colour <- c("blue", "green", "red")
colnum <- rep(1:3, n/3+1)[1:n] - 1

data <- data.frame(x=x, colour=colour[colnum+1])

glmformula <- formula( ~ x + I(x^2) + colour)
glmfamily <- binomial(link="logit")
Xdesign0 <- model.matrix(glmformula, data=data)

#    multiple subjects
Xdesign <- NULL
for (i in 1:N) Xdesign <- rbind(Xdesign, Xdesign0)

# --- Parameter Values and Simulation ---

Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

delta <- c(0.5, 0.5)

beta <- matrix(c(-1, -1.2,
                 -2, -1.8,
                  3,  2.8,
                  1,  0.8, 
                  2,  2.2), 
               ncol=ncol(Pi), nrow=ncol(Xdesign), byrow=TRUE)

y <- mmglmlong1(NULL, Pi, delta, glmfamily, beta, Xdesign, sigma=sd,
                longitude=rep(1:N, each=n), size=rep(200, N*n))

y <- simulate(y, seed=5)

# --- Estimation ---

#   Increase maxiter below to achieve convergence
#   Has been restricted to minimise time of package checks
tmp <- BaumWelch(y, bwcontrol(posdiff=FALSE, maxiter=1))
print(summary(tmp))
}

\keyword{classes}

