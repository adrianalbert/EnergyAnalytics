\name{probhmm}
\alias{probhmm}

\title{Conditional Distribution Function of DTHMM}

\description{
Calculates the distribution function at each point for a \code{\link{dthmm}} process given the complete observed process except the given point.
}

\usage{
probhmm(logalpha, logbeta, Pi, delta, cumprob)
}

\arguments{
  \item{logalpha}{an \eqn{n \times m}{n*m} matrix containing the logarithm of the forward probabilities.}
  \item{logbeta}{an \eqn{n \times m}{n*m} matrix containing the logarithm of the backward probabilities.}
  \item{Pi}{is the \eqn{m \times m}{m*m} transition probability matrix of the hidden Markov chain.}
  \item{delta}{is the marginal probability distribution of the \eqn{m} hidden states at the first time point.}
  \item{cumprob}{an \eqn{n \times m}{n*m} matrix where the \eqn{(i,k)}th element is \eqn{\Pr\{ X_i \le x_i \,|\, C_k = c_k\}}{Pr{Xi <= xi | Ck = ck}}.}
}

\details{
Let \eqn{X^{(-i)}} denote the entire process, except with the point \eqn{X_i}{Xi} removed. The distribution function at the point \eqn{X_i}{Xi} is
\deqn{
\Pr\{ X_i \le x_i \,|\, X^{(-i)} = x^{(-i)} \}\,.
}{
Pr{ Xi <= xi | X^{(-i)} = x^{(-i)} } .
}
This \R function calculates the distribution function for each point \eqn{X_i}{Xi} for \eqn{i=1, \cdots, n}{i=1, ..., n}. This is done by using the forward and backward probabilities before and after the \eqn{i}th point, respectively.

In the programming code, note the subtraction of the mean. This is to stop underflow when the exponential is taken. Removal of the mean is automatically compensated for by the fact that the same factor is removed in both the numerator and denominator.
}

\value{
A vector containing the probability.
}


\seealso{\code{\link{residuals}}}

\keyword{distribution}
