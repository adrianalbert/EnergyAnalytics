\name{forwardback}
\alias{forwardback}
\alias{forward}
\alias{backward}
\alias{forwardback.dthmm}

\title{Forward and Backward Probabilities of DTHMM}

\description{
These functions calculate the forward and backward probabilities for a \code{\link{dthmm}} process, as defined in MacDonald & Zucchini (1997, Page 60).
}

\usage{
backward(x, Pi, distn, pm, pn = NULL)
forward(x, Pi, delta, distn, pm, pn = NULL)
forwardback(x, Pi, delta, distn, pm, pn = NULL, fortran = TRUE)
forwardback.dthmm(Pi, delta, prob, fortran = TRUE, fwd.only = FALSE)
}

\arguments{
  \item{x}{is a vector of length \eqn{n} containing the observed process.}
  \item{Pi}{is the \eqn{m \times m}{m*m} transition probability matrix of the hidden Markov chain.}
  \item{delta}{is the marginal probability distribution of the \eqn{m} hidden states.}
  \item{distn}{is a character string with the distribution name, e.g. \code{"norm"} or \code{"pois"}. If the distribution is specified as \code{"wxyz"} then a probability (or density) function called \code{"dwxyz"} should be available, in the standard \R format (e.g. \code{\link{dnorm}} or \code{\link{dpois}}).}
  \item{pm}{is a list object containing the current (Markov dependent) parameter estimates associated with the distribution of the observed process (see \code{\link{dthmm}}).}
  \item{pn}{is a list object containing the observation dependent parameter values associated with the distribution of the observed process (see \code{\link{dthmm}}).}
  \item{prob}{an \eqn{n \times m}{n*m} matrix containing the observation probabilities or densities (rows) by Markov state (columns).}
  \item{fortran}{logical, if \code{TRUE} (default) use the Fortran code, else use the \R code.}
  \item{fwd.only}{logical, if \code{FALSE} (default) calculate both forward and backward probabilities; else calculate and return only forward probabilities and log-likelihood.}
}

\value{
The function \code{forwardback} returns a list with two matrices containing the forward and backward (log) probabilities, \code{logalpha} and \code{logbeta}, respectively, and the log-likelihood (\code{LL}).

The functions \code{backward} and \code{forward} return a matrix containing the forward and backward (log) probabilities, \code{logalpha} and \code{logbeta}, respectively.
}

\details{
Denote the \eqn{n \times m}{n*m} matrices containing the forward and backward probabilities as \eqn{A} and \eqn{B}, respectively. Then the \eqn{(i,j)}th elements are
\deqn{
\alpha_{ij} = \Pr\{ X_1 = x_1, \cdots, X_i = x_i, C_i = j \}
}{
alpha_{ij} = Pr{ X_1 = x_1, ..., X_i = x_i, C_i = j }
}
and
\deqn{
\beta_{ij} = \Pr\{ X_{i+1} = x_{i+1}, \cdots, X_n = x_n \,|\, C_i = j \} \,.
}{
beta_{ij} = Pr{ X_{i+1} = x_{i+1}, ..., X_n = x_n | C_i = j } .
}
Further, the diagonal elements of the product matrix \eqn{A B^\prime}{AB'} are all the same, taking the value of the log-likelihood.
}

\author{The algorithm has been taken from Zucchini (2005).}

\seealso{
\code{\link{logLik}}
}

\references{
Cited references are listed on the \link{HiddenMarkov} manual page.
}

\examples{
#    Set Parameter Values

Pi <- matrix(c(1/2, 1/2,   0,   0,   0,
               1/3, 1/3, 1/3,   0,   0,
                 0, 1/3, 1/3, 1/3,   0,
                 0,   0, 1/3, 1/3, 1/3,
                 0,   0,   0, 1/2, 1/2),
             byrow=TRUE, nrow=5)

p <- c(1, 4, 2, 5, 3)
delta <- c(0, 1, 0, 0, 0)

#------   Poisson HMM   ------

x <- dthmm(NULL, Pi, delta, "pois", list(lambda=p), discrete=TRUE)

x <- simulate(x, nsim=10)

y <- forwardback(x$x, Pi, delta, "pois", list(lambda=p))

# below should be same as LL for all time points
print(log(diag(exp(y$logalpha) \%*\% t(exp(y$logbeta)))))
print(y$LL)

#------   Gaussian HMM   ------

x <- dthmm(NULL, Pi, delta, "norm", list(mean=p, sd=p/3))

x <- simulate(x, nsim=10)

y <- forwardback(x$x, Pi, delta, "norm", list(mean=p, sd=p/3))

# below should be same as LL for all time points
print(log(diag(exp(y$logalpha) \%*\% t(exp(y$logbeta)))))
print(y$LL)
}

\keyword{distribution}
