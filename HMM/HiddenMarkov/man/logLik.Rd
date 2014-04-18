\name{logLik}
\alias{logLik}
\alias{logLik.dthmm}
\alias{logLik.mmglm0}
\alias{logLik.mmglm1}
\alias{logLik.mmglmlong1}
\alias{logLik.mmpp}

\title{Log Likelihood of Hidden Markov Model}
\description{
Provides methods for the generic function \code{\link[stats]{logLik}}.
}

\usage{
\method{logLik}{dthmm}(object, fortran=TRUE, ...)
\method{logLik}{mmglm0}(object, fortran=TRUE, ...)
\method{logLik}{mmglm1}(object, fortran=TRUE, ...)
\method{logLik}{mmglmlong1}(object, fortran=TRUE, ...)
\method{logLik}{mmpp}(object, fortran=TRUE, ...)
}

\arguments{
  \item{object}{an object with class \code{"\link{dthmm}"}, \code{"\link{mmglm0}"}, \code{"\link{mmglm1}"}, \code{"\link{mmglmlong1}"} or \code{"\link{mmpp}"}.}
  \item{fortran}{logical, if \code{TRUE} (default) use the Fortran code, else use the \R code.}
  \item{...}{other arguments.}
}

\details{
The methods provided here will always recalculate the log-likelihood even if it is already contained within the \code{object}. This enables the user to change parameter or data values within the \code{object} and recalculate the log-likelihood for the revised configuration.

The code for the methods \code{"\link{dthmm}"}, \code{"\link{mmglm0}"}, \code{"\link{mmglm1}"}, \code{"\link{mmglmlong1}"} and \code{"\link{mmpp}"} can be viewed by typing \code{logLik.dthmm}, \code{logLik.mmglm0}, \code{logLik.mmglm1}, \code{logLik.mmglmlong1} or \code{logLik.mmpp}, respectively, on the \R command line.
}

\value{
Returns the value of the log-likelihood.
}

\examples{
Pi <- matrix(c(1/2, 1/2,   0,
               1/3, 1/3, 1/3,
                 0, 1/2, 1/2),
             byrow=TRUE, nrow=3)

x <- dthmm(NULL, Pi, c(0,1,0), "norm",
           list(mean=c(1, 6, 3), sd=c(1, 0.5, 1)))

x <- simulate(x, nsim=100)

print(logLik(x))
}

\keyword{methods}
