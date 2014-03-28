\name{viterbi}
\alias{viterbi}
\title{Viterbi algorithm}
\description{This function calculates the optimal hidden states sequence using Viterbi's algorithm.}
\usage{
viterbi(HMM, obs)
}
\arguments{
    \item{HMM}{a HMMClass or a HMMFitClass object.}
    \item{obs}{The vector, matrix, data frame, list of vectors or list of matrices of observations.}
}
\value{a viterbiClass object which is a list with:
\item{States}{Sequence of hidden states in 1...nStates}
\item{logViterbiScore}{logarithm of the Viterbi's Score.}
\item{logProbSeq}{logarithm of probability of having the sequence of states conditionally to having the observations.}
}
\examples{
 data(n1d_3s)
 ResFit <- HMMFit(obs_n1d_3s, nStates=3)
 VitPath <- viterbi(ResFit, obs_n1d_3s)
}
\references{
    Among hundreds of tutorials, you can have a look to use \cr
    Phil Blunsom (2004) \emph{ Hidden Markov Models. }
    \url{http://www.cs.mu.oz.au/460/2004/materials/hmm-tutorial.pdf}}

\seealso{\code{\link{HMMSet}}, \code{\link{HMMFit}}}
\keyword{hplot}
