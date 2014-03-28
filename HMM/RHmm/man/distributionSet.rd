\name{distributionSet}
\alias{distributionSet}
\alias{distributionClass}
\title{Set the parameters for the distributions of observations}
\description{This function is used to create a distributionClass object which contains the parameters of the distribution of the observations
    for each hidden state. Since distributions can be univariate or multivariate, discrete or continuous,
    the different values of a distributionClass object depend of the nature of the distribution.}


\usage{
distributionSet(dis, ...)
}

\arguments{
    \item{dis}{Name of the distribution of observations. In 'NORMAL', 'DISCRETE', 'MIXTURE'.}
    \item{...}{Other parameters. See details.}
}

\details{
Typical usages are:
\itemize{
    \item distributionSet(dis="NORMAL", mean, var)
    \item distributionSet(dis="NORMAL", mean, cov)
    \item distributionSet(dis="MIXTURE", mean, var, proportion)
    \item distributionSet(dis="MIXTURE", mean, cov, proportion)
    \item distributionSet(dis="DISCRETE", proba, labels=NULL)
}
The parameters are:
\describe{
    \item{mean}{
        \emph{- Univariate normal distribution: }a vector of the means for each state of the model.\cr
        \emph{- Multivariate normal distribution: }a list of the mean vectors for each state of the model.\cr
        \emph{- Mixture of univariate normal distribution: }a list of vectors of the mixture means for each state of the model.\cr
        \emph{- Mixture of multivariate normal distribution: }a list of lists of vectors of means for each state and each component of the mixture of the model.
        }
    \item{var}{
        \emph{- Univariate normal distribution: }a vector of the variances for each states of the model.\cr
        \emph{- Mixture of univariate normal distribution: }a list of vectors of the mixture variances for each states of the model.
        }
    \item{cov}{
        \emph{- Multivariate normal distribution: }a list of covariance matrices of the multivariate normal distribution for each state of the model\cr
        \emph{- Mixture of multivariate normal distribution: }a list of list of covariance matrices for each state and each component of the mixture.
         }
    \item{proportion}{A list of vector of the mixture proportions for each state of the model.}
    \item{proba}{A list of vector of discrete probabilities for each state of the model.}
    \item{labels}{A vector of the labels of the discrete observations. Default NULL.}
    }
}

\value{An \sQuote{distributionClass} class object with some of the following elements:
    \describe{
        \item{dis}{The name of the distribution.}
        \item{nStates}{Number of hidden states.}
        \item{dimObs}{Dimension of observations.}
        \item{nMixt}{Number of mixtures for mixture of normal distributions.}
        \item{nLevels}{Number of levels for discrete distributions.}
        \item{mean}{The \sQuote{mean} argument for univariate normal, mixture of univariate normal and multivariate normal distributions.}
        \item{var}{The \sQuote{var} argument for univariate normal and mixture of univariate normal distributions}
        \item{cov}{The \sQuote{cov} argument for multivariate normal and mixture of multivaiate normal distributions}
        \item{proba}{The \sQuote{proba} argument for discrete distributions}
        }
    }


\examples{
    # 3 hidden states Markov Model with univariate normal distributions
    # for the observations
    #   obs | hidden state = 1 are N(1, 1)
    #   obs | hidden state = 2 are N(-2, 2)
    #   obs | hidden state = 3 are N(5, 4)
        n_1d_3s <- distributionSet("NORMAL", mean=c(1, -2, 5), var=c(1, 2, 4))
    # 2 hidden states Markov Model with bivariate normal distributions
    # for the observations
    #   obs | hidden state = 1 are N(m1, cov1)
    #   obs | hidden state = 2 are N(m2, cov2)
        m1 <- c(1,1)
        m2 <- c(-2, -2)
        cov1 <- matrix(c(1, 1, 1, 4), nrow=2)
        cov2 <- matrix(c(1, -1, -1, 9), nrow=2)
        n_2d_2s <- distributionSet("NORMAL", mean=list(m1, m2),
                                        cov=list(cov1, cov2))
    # 3 hidden states Markov Model with a mixture of two normal
    # distributions for the observations
    # obs | hidden state = i are:
    #   pi[1] * N(mmi[1], vari[1]) + pi[2] * N(mmi[2], vari[2])

        mm1 <- c(1, -1)
        mm2 <- c(-2, 2)
        mm3 <- c(5, 5)
        var1 <- c(1, 2)
        var2 <- c(2, 3)
        var3 <- c(1, 1)
        p1 <- c(0.5, 0.5)
        p2 <- c(0.8, 0.2)
        p3 <- c(0.3, 0.7)
        mn_2s <- distributionSet("MIXTURE", mean=list(mm1, mm2, mm3),
            var=list(var1, var2, var3), proportion=list(p1, p2, p3))
    # 2 hidden states Markov Model with discrete observations
        dp1 <- c(0.2, 0.3, 0.3, 0.2)
        dp2 <- c(0.1, 0.1, 0.1, 0.7)
        labels <- c("I", "M", "A", "G")
        d_2s <- distributionSet("DISCRETE", proba=list(dp1, dp2),
                                labels=labels)
    # 2 hidden states Markov model with mixture of 3 2-d gaussian distribution
        q1 <- rep(1/3, 3)
        q2 <- runif(3)
        q2 <- q2/sum(q2)
        cov3 <- matrix(c(1,2,2,10), nrow=2)
        cov4 <- matrix(c(1, 0, 0, 1), nrow=2)
        cov5 <- matrix(c(2,4,4,50), nrow=2)
        cov6 <- matrix(c(25,1, 1, 2), nrow=2)
        mm4 <- c(100, 20)
        mm5 <- c(20, -20)
        mm6 <- c(0, 0)
        m_2d_2s <- distributionSet("MIXTURE", mean=list(list(mm1,mm2,mm3), list(mm4,mm5,mm6)),
            cov=list(list(cov1,cov2,cov3), list(cov4,cov5,cov6)), proportion=list(q1,q2))

}
\keyword{models}
\keyword{distribution}
