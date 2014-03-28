\name{asymptoticIterSimCovMat
}
\alias{asymptoticIterSimCovMat}
\title{Compute the asymptotic covariance matrix of a fitted HMM by simulation}
\description{This \sQuote{new} function computes the empirical asymptotic covariance matrix of the fitted HMM.}
\usage{
asymptoticIterSimCovMat(HMM, obs, nSimul, verbose=FALSE, oldCovMat=NULL)
}
\arguments{
    \item{HMM}{a HMMClass or HMMFitClass object}
    \item{obs}{A vector, a matrix, a data frame, a list of vectors or a list of matrices of observations. \bold{See HMMFit}.}
    \item{nSimul}{The number of simulation}
    \item{verbose}{A boolean. if true, displays some informations. Default false.}
    \item{oldCovMat}{An object containing \itemize{
        \item{esp2man: }{The current matrix of the empirical mean of \eqn{\theta\%*\%t(\theta)}}
        \item{mean: }{The current vector of the empirical mean of \eqn{\theta}}
        \item{cov: }{the current empirical covariance matrix of \eqn{\theta}}
        \item{nSimul: }{The current number of simulations}
        }
        where \eqn{\theta} is the vector of all parameters of the HMM.
    }
}
\value{An object with the same attributes than \sQuote{oldCovMat} parameter.}


\section{Numerical computations}{This is an ``experimental'' method. The HMM model is simulated nSimul times then fitted
and the empirical covariance matrix is computed.}

\examples{
    # Fit a 3 states 1D-gaussian model
    data(n1d_3s)
    Res <- HMMFit(obs_n1d_3s, nStates=3)
    # First 10 computations of covariance matrix
    Cov <- asymptoticIterSimCovMat(Res, obs_n1d_3s, 10)
    # 10 more computations of covariance matrix
    Cov <- asymptoticIterSimCovMat(Res, obs_n1d_3s, 10, verbose=TRUE, oldCovMat=Cov)
    Res<-setAsymptoticCovMat(Res, Cov$cov)
    summary(Res)
    }

\seealso{\code{\link{setAsymptoticCovMat}}, \code{\link{asymptoticCov}}}
