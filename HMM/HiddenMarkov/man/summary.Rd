\name{summary}
\alias{summary}
\alias{summary.dthmm}
\alias{summary.mmpp}
\alias{summary.mmglm0}
\alias{summary.mmglm1}
\alias{summary.mmglmlong1}

\title{Summary of Hidden Markov Model}

\description{
Provides methods for the generic function \code{\link[base]{summary}}.
}

\usage{
\method{summary}{dthmm}(object, ...)
\method{summary}{mmglm0}(object, ...)
\method{summary}{mmglm1}(object, ...)
\method{summary}{mmglmlong1}(object, ...)
\method{summary}{mmpp}(object, ...)
}

\arguments{
  \item{object}{an object with class \code{"\link{dthmm}"}, \code{"\link{mmglm0}"}, \code{"\link{mmglm1}"}, \code{"\link{mmglmlong1}"} or \code{"\link{mmpp}"}.}
  \item{...}{other arguments.}
}


\details{
The code for the methods \code{"\link{dthmm}"}, \code{"\link{mmglm0}"}, \code{"\link{mmglm1}"}, \code{"\link{mmglmlong1}"} and \code{"\link{mmpp}"} can be viewed by typing \code{summary.dthmm}, \code{summary.mmglm0}, \code{summary.mmglm1},\cr \code{summary.mmglmlong1} or \code{summary.mmpp}, respectively, on the \R command line.
}

\value{
A list object with a reduced number of components, mainly the parameter values.
}

\examples{
Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

x <- dthmm(NULL, Pi, c(0, 1), "beta",
           list(shape1=c(2, 6), shape2=c(6, 2)))

x <- simulate(x, nsim=2000)

print(summary(x))
}

\keyword{methods}

