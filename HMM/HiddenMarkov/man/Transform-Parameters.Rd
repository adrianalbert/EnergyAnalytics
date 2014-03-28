\name{Transform-Parameters}
\alias{Transform-Parameters}
\alias{Pi2vector}
\alias{vector2Pi}
\alias{Q2vector}
\alias{vector2Q}

\title{Transform Transition or Rate Matrices to Vector}

\description{
These functions transform \eqn{m \times m}{m*m} transition probability matrices or \eqn{Q} matrices to a vector of length \eqn{m(m-1)}, and back. See Details.
}

\usage{
Pi2vector(Pi)
vector2Pi(p)

Q2vector(Q)
vector2Q(p)
}

\arguments{
  \item{Pi}{an \eqn{m} by \eqn{m} transition probability matrix.}
  \item{Q}{an \eqn{m} by \eqn{m} rate matrix.}
  \item{p}{a vector of length \eqn{m(m-1)}.}
}

\details{
The function \code{Pi2vector} maps the \eqn{m} by \eqn{m} transition probability matrix of a discrete time HMM to a vector of length \eqn{m(m-1)}, and \code{vector2Pi} has the reverse effect. They use a logit like transformation so that the parameter space is on the whole real line thus avoiding hard boundaries which cause problems for many optimisation procedures (see \code{\link{neglogLik}}).

Similarly, the function \code{Q2vector} maps the \eqn{m} by \eqn{m} rate matrix \eqn{Q} of an MMPP process to a vector of length \eqn{m(m-1)}, and \code{vector2Q} has the reverse effect. They use a log transformation so that the parameter space is on the whole real line thus avoiding hard boundaries which cause problems for many optimisation procedures (see \code{\link{neglogLik}}).
}

\value{
The functions \code{Pi2vector} and \code{Q2vector} return a vector of length \eqn{m(m-1)}, the function \code{vector2Pi} returns an \eqn{m} by \eqn{m} transition probability matrix, and \code{vector2Q} returns an \eqn{m \times m}{m*m} rate matrix \eqn{Q}.
}

\seealso{\code{\link{neglogLik}}}

\examples{
Pi <- matrix(c(0.8, 0.1, 0.1,
               0.1, 0.6, 0.3,
               0.2, 0.3, 0.5),
             byrow=TRUE, nrow=3)

print(vector2Pi(Pi2vector(Pi)))

#------------------------------------------------

Q <- matrix(c(-8,  5,  3,
               1, -4,  3,
               2,  5, -7),
            byrow=TRUE, nrow=3)

print(vector2Q(Q2vector(Q)))
}

\keyword{optimize}
\keyword{misc}
