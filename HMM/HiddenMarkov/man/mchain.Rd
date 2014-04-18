\name{mchain}
\alias{mchain}

\title{Markov Chain Object}

\description{
Creates a Markov chain object with class \code{"mchain"}. It does not simulate data.
}

\usage{
mchain(x, Pi, delta, nonstat = TRUE)
}

\arguments{
  \item{x}{is a vector of length \eqn{n} containing the observed process, else it is specified as \code{NULL}. This is used when there are no data and a process is to be simulated.}
  \item{Pi}{is the \eqn{m \times m}{m*m} transition probability matrix of the Markov chain.}
  \item{delta}{is the marginal probability distribution of the \eqn{m} state Markov chain at the first time point.}
  \item{nonstat}{is logical, \code{TRUE} if the homogeneous Markov chain is assumed to be non-stationary, default. See \dQuote{Details} below.}
}

\value{
A \code{\link[base]{list}} object with class \code{"mchain"}, containing the above arguments as named components.
}

\examples{
Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

#    Create a Markov chain object with no data (NULL)
x <- mchain(NULL, Pi, c(0,1))

#    Simulate some data
x <- simulate(x, nsim=2000)

#   estimate transition probabilities
estPi <- table(x$mc[-length(x$mc)], x$mc[-1])
rowtotal <- estPi \%*\% matrix(1, nrow=nrow(Pi), ncol=1)
estPi <- diag(as.vector(1/rowtotal)) \%*\% estPi
print(estPi)
}


\keyword{classes}
