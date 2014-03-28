\name{Estep}
\alias{Estep}

\title{E-Step of EM Algorithm for DTHMM}

\description{
Performs the \emph{expectation} step of the EM algorithm for a \code{\link{dthmm}} process. This function is called by the \code{\link{BaumWelch}} function. The Baum-Welch algorithm referred to in the HMM literature is a version of the EM algorithm.
}

\usage{
Estep(x, Pi, delta, distn, pm, pn = NULL)
}

\arguments{
  \item{x}{is a vector of length \eqn{n} containing the observed process.}
  \item{Pi}{is the current estimate of the \eqn{m \times m}{m*m} transition probability matrix of the hidden Markov chain.}
  \item{distn}{is a character string with the distribution name, e.g. \code{"norm"} or \code{"pois"}. If the distribution is specified as \code{"wxyz"} then a probability (or density) function called \code{"dwxyz"} should be available, in the standard \R format (e.g. \code{\link{dnorm}} or \code{\link{dpois}}).}
  \item{pm}{is a list object containing the current (Markov dependent) parameter estimates associated with the distribution of the observed process (see \code{\link{dthmm}}).}
  \item{pn}{is a list object containing the observation dependent parameter values associated with the distribution of the observed process (see \code{\link{dthmm}}).}
  \item{delta}{is the current estimate of the marginal probability distribution of the \eqn{m} hidden states.}
}

\value{
A \code{list} object is returned with the following components.
  \item{u}{an \eqn{n \times m}{n*m} matrix containing estimates of the conditional expectations. See \dQuote{Details}.}
  \item{v}{an \eqn{n \times m \times m}{n*m*m} array containing estimates of the conditional expectations. See \dQuote{Details}.}
  \item{LL}{the current value of the log-likelihood.}
}

\details{
Let \eqn{u_{ij}} be one if \eqn{C_i=j} and zero otherwise. Further, let \eqn{v_{ijk}} be one if \eqn{C_{i-1}=j} and \eqn{C_i=k}, and zero otherwise. Let \eqn{X^{(n)}} contain the complete observed process. Then, given the current model parameter estimates, the returned value \code{u[i,j]} is
\deqn{
\widehat{u}_{ij} = \mbox{E}[u_{ij} \, | \, X^{(n)}] = \Pr\{C_i=j \, | \, X^{(n)} = x^{(n)} \} \,,
}{
hat{u}_{ij} = E[u_{ij} | X^{(n)}] = Pr{ C_i=j | X^{(n)} = x^{(n)} } ,
}
and \code{v[i,j,k]} is
\deqn{
\widehat{v}_{ijk} = \mbox{E}[v_{ijk} \, | \, X^{(n)}] = \Pr\{C_{i-1}=j, C_i=k \, | \, X^{(n)} = x^{(n)} \}\,,
}{
hat{v}_{ijk} = E[v_{ijk} | X^{(n)}] = Pr{ C_{i-1}=j, C_i=k | X^{(n)} = x^{(n)} } ,
}
where \eqn{j,k = 1, \cdots, m}{j,k = 1, ..., m} and \eqn{i = 1, \cdots, n}{i = 1, ..., n}.
}


\author{The algorithm has been taken from Zucchini (2005).}

\seealso{
\code{\link{BaumWelch}}, \code{\link{Mstep}}
}

\references{
Cited references are listed on the \link{HiddenMarkov} manual page.
}

\keyword{optimize}

