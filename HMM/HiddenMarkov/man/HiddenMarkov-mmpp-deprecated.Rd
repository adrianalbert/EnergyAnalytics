\name{HiddenMarkov-mmpp-deprecated}
\alias{HiddenMarkov-mmpp-deprecated}
\alias{backward0.mmpp}
\alias{forward0.mmpp}
\alias{logLikmmpp}
\alias{Estep0.mmpp}
\alias{Baum.Welch.mmpp}
\alias{Baum.Welch0.mmpp}
\alias{sim.mmpp}

\title{Markov Modulated Poisson Process - Deprecated Functions}

\description{
These functions are deprecated and will ultimately be removed from the package. Please change to the revised versions: \code{\link{BaumWelch}}, \code{\link{Estep.mmpp}}, \code{\link{forwardback.mmpp}}, \code{\link[HiddenMarkov]{simulate}} or \code{\link{logLik}}.
}

\usage{
backward0.mmpp(tau, Q, lambda)
forward0.mmpp(tau, Q, delta, lambda)

logLikmmpp(tau, Q, delta, lambda)

Estep0.mmpp(tau, Q, delta, lambda)

Baum.Welch.mmpp(tau, Q, delta, lambda, nonstat = TRUE,
                maxiter = 500, tol = 1e-05, prt = TRUE,
                converge = expression(diff < tol))
Baum.Welch0.mmpp(tau, Q, delta, lambda, nonstat = TRUE,
                 maxiter = 500, tol = 1e-05, prt = TRUE,
                 converge = expression(diff < tol))

sim.mmpp(n, initial, Q, lambda)
}

\arguments{
  \item{tau}{vector containing the interevent times. Note that the first event is at time zero.}
  \item{Q}{the infinitesimal generator matrix of the Markov process.}
  \item{lambda}{a vector containing the Poisson rates.}
  \item{delta}{is the marginal probability distribution of the \eqn{m} hidden states at time zero.}
  \item{n}{number of Poisson events to be simulated.}
  \item{initial}{integer, being the initial hidden Markov state \eqn{(1, \cdots, m)}.}
  \item{nonstat}{is logical, \code{TRUE} if the homogeneous Markov chain is assumed to be non-stationary, default.}
  \item{maxiter}{is the maximum number of iterations, default is 500.}
  \item{tol}{is the convergence criterion, being the difference between successive values of the log-likelihood; default is 0.00001.}
  \item{prt}{is logical, and determines whether information is printed at each iteration; default is \code{TRUE}.}
  \item{converge}{is an expression giving the convergence criterion.}
}

\details{
The functions with a suffix of zero are non-scaled, and hence will have numerical problems for series containing larger numbers of events; and are \emph{much} slower.

These functions use the algorithm given by Ryden (1996) based on eigenvalue decompositions.
}

\keyword{internal}

