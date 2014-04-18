\name{compdelta}
\alias{compdelta}

\title{Marginal Distribution of Stationary Markov Chain}

\description{
Computes the marginal distribution of a \emph{stationary} Markov chain with transition probability matrix \eqn{\Pi}{Pi}. The \eqn{m} discrete states of the Markov chain are denoted by \eqn{1, \cdots, m}{1, ..., m}.
}

\usage{
compdelta(Pi)
}

\arguments{
  \item{Pi}{is the \eqn{m \times m}{m*m} transition probability matrix of the Markov chain.}
}

\details{
If the Markov chain is stationary, then the marginal distribution \eqn{\delta}{delta} satisfies
\deqn{
\delta = \delta \Pi \,.
}{
delta = delta Pi.
}
Obviously,
\deqn{
\sum_j^m \delta_j = 1.
}{
sum_j^m delta_j = 1.
}
}

\value{
A numeric vector of length \eqn{m} containing the marginal probabilities.
}

\examples{
Pi <- matrix(c(1/2, 1/2,   0,   0,   0,
               1/3, 1/3, 1/3,   0,   0,
                 0, 1/3, 1/3, 1/3,   0,
                 0,   0, 1/3, 1/3, 1/3,
                 0,   0,   0, 1/2, 1/2),
             byrow=TRUE, nrow=5)

print(compdelta(Pi))
}
\keyword{distribution}
