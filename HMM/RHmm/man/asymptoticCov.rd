\name{asymptoticCov}
\alias{asymptoticCov}
\title{Asymptotic covariance matrix of the HMM parameters}
\description{This function calculates the empirical asymptotic covariance matrix of the HMM parameters}
\usage{
asymptoticCov(HMM, obs)
}
\arguments{
    \item{HMM}{A HMMClass or a HMMFitClass object}
    \item{obs}{The vector, matrix, data frame, list of vectors or list of matrices of observations}
}
\value{A matrix}

\section{Numerical computations}{
   The Information matrix (of the independent parameters) is computed using the Lystig and Hugues's algorithm. Then the covariance matrix is computed by inversion of this information matrix.
}

\examples{
  data(n1d_3s)
  Res_n1d_3s<-HMMFit(obs_n1d_3s, nStates=3)
  covMat <- asymptoticCov(Res_n1d_3s, obs_n1d_3s)
}

\references{
     Lystig Theodore C. and Hugues James P. (2002) \emph{Exact Computation of the Observed Information Matrix for Hidden Markov Models}, Journal of Computational and Graphical Statistics, Vol. 11, No 3, 678-689.

}

\seealso{HMMFit}
