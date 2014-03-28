\name{computeScoreAndInformation}
\alias{computeScoreAndInformation}
\title{Score and Information matrix of the HMM parameters}
\description{This function calculates the score and the information matrix of the independent parameters of the HMM, using Lystig and Hugues's algorithm.}
\usage{
computeScoreAndInformation(HMM, obs) }
\arguments{
    \item{HMM}{A HMMClass or a HMMFitClass object}
    \item{obs}{The vector, matrix, data frame, list of vectors or list of matrices of observations}
}

\value{
    \item{score}{the score vector of the independent parameters}
    \item{information}{the information matrix of the independent parameters}
}


\examples{
  data(n1d_3s)
  Res_n1d_3s<-HMMFit(obs_n1d_3s, nStates=3)
  ScoreAndInfor <- computeScoreAndInformation(Res_n1d_3s, obs_n1d_3s)
}

\references{
     Lystig Theodore C. and Hugues James P. (2002) \emph{Exact Computation of the Observed Information Matrix
     for Hidden Markov Models}, Journal of Computational and Graphical Statistics, Vol. 11, No 3, 678-689.

}

\seealso{HMMFit}
