\name{HiddenMarkov-dthmm-deprecated}
\alias{HiddenMarkov-dthmm-deprecated}
\alias{Baum.Welch}
\alias{sim.hmm}
\alias{sim.hmm1}
\alias{sim.markov}
\alias{Viterbihmm}

\title{Discrete Time HMM - Deprecated Functions}

\description{
These functions are deprecated and will ultimately be removed from the package. Please change to the object orientated versions: \code{\link{BaumWelch}}, \code{\link[HiddenMarkov]{residuals}}, \code{\link[HiddenMarkov]{simulate}} or \code{\link{Viterbi}}.
}

\usage{
Baum.Welch(x, Pi, delta, distn, pm, pn = NULL, nonstat = TRUE,
           maxiter = 500, tol = 1e-05, prt = TRUE,
           posdiff = (distn[1]!="glm"))
residualshmm(x, Pi, delta, distn, pm, pn = NULL, discrete = FALSE)
sim.hmm(n, initial, Pi, distn, pm, pn = NULL)
sim.hmm1(n, initial, Pi, distn, pm)
sim.markov(n, initial, Pi)
Viterbihmm(x, Pi, delta, distn, pm, pn = NULL)
}

\arguments{
  \item{x}{is a vector of length \eqn{n} containing the observed process.}
  \item{n}{length of process.}
  \item{initial}{integer, being the initial hidden Markov state \eqn{(1, \cdots, m)}.}
  \item{Pi}{is the \eqn{m \times m}{m*m} transition probability matrix of the hidden Markov chain.}
  \item{delta}{is the marginal probability distribution of the \eqn{m} hidden states at the first time point.}
  \item{distn}{is a character string with the distribution name, e.g. \code{"norm"} or \code{"pois"}. If the distribution is specified as \code{"wxyz"} then a distribution function called \code{"pwxyz"} should be available, in the standard \R format (e.g. \code{\link{pnorm}} or \code{\link{ppois}}).}
  \item{pm}{is a list object containing the (Markov dependent) parameter values associated with the distribution of the observed process (see \code{\link{dthmm}}).}
  \item{pn}{is a list object containing the observation dependent parameter values associated with the distribution of the observed process (see \code{\link{dthmm}}).}
  \item{discrete}{is logical, and is \code{TRUE} if \code{distn} is a discrete distribution.}
  \item{nonstat}{is logical, \code{TRUE} if the homogeneous Markov chain is assumed to be non-stationary, default. See \dQuote{Details} below.}
  \item{maxiter}{is the maximum number of iterations, default is 500.}
  \item{tol}{is the convergence criterion, being the difference between successive values of the log-likelihood; default is 0.00001.}
  \item{prt}{is logical, and determines whether information is printed at each iteration; default is \code{TRUE}.}
  \item{posdiff}{is logical, and determines whether the iterative process stops if a negative log-likelihood difference occurs.}
}

\details{
The function \code{sim.hmm1} will run faster for cases where the argument \code{pn} is \code{NULL}.
}

\keyword{internal}
