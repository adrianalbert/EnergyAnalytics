\name{HMMSet}
\alias{HMMSet}
\alias{HMMClass}
\title{Set the parameters for the hidden Markov models}

\description{This function is used to create a HMMClass object which contains the parameters of the HMM. An HMM is described by an initial state probability vector,
transition matrices and a distributionClass object. }

\usage{
    HMMSet(initProb, transMat, ...)
}
\arguments{
    \item{initProb}{The vector of probabilities of the initial state}
    \item{transMat}{The transition matrix of the hidden Markov chain. Since version 1.3.4 of RHmm, a list of transition matrices can be specified
in order to define an inhomogeneous HMM.}
    \item{...}{Other parameters. See details.}
}

\details{
Typical usages are:\cr
\itemize{
    \item HMMSet(initProb, transMat, distribution)
    \item HMMSet(initProb, transMat, dis="NORMAL", mean, var)
    \item HMMSet(initProb, transMat, dis="NORMAL", mean, cov)
    \item HMMSet(initProb, transMat, dis="MIXTURE", mean, var, proportion)
    \item HMMSet(initProb, transMat, dis="DISCRETE", proba, labels=NULL)
}
The different arguments are:
    
\describe{
    \item{distribution}{The distributionClass object of the observations}
    \item{dis}{dis parameter. See \bold{distributionSet}}
    \item{mean}{mean parameter. See \bold{distributionSet}}
    \item{var}{var parameter. See \bold{distributionSet}}
    \item{cov}{cov parameter. See \bold{distributionSet}}
    \item{proportion}{proportion parameter. See \bold{distributionSet}}
    \item{proba}{proba parameter. See \bold{distributionSet}}
    \item{labels}{labels parameter. See \bold{distributionSet}}
    }
}

\value{An \sQuote{HMMClass} class object with the following elements:
    \describe{
        \item{initProb}{Initial state probabilities vector}
        \item{transMat}{Transition matrix}
        \item{distribution}{The distributionClass object which describes the conditional distribution. See \bold{distributionSet.}}
    }
}
\examples{
    # 3 hidden states Markov Model with univariate normal distributions
    # for the observations
    #   obs | hidden state = 1 are N(1, 1)
    #   obs | hidden state = 2 are N(-2, 2)
    #   obs | hidden state = 3 are N(5, 4)

        n_1d_3s <- distributionSet("NORMAL", c(1, -2, 5), c(1, 2, 4))
        initProb3 <- rep(1,3)/3
        transMat3 <- rbind(c(0.5, 0.4, 0.1), c(0.3, 0.4, 0.3),
            c(0.2, 0.1, 0.7))
        hmm1 <- HMMSet(initProb3, transMat3, n_1d_3s)
        # or directly
        hmm2 <- HMMSet(initProb3, transMat3, "NORMAL", mean=c(1, -2, 5),
            var=c(1, 2, 4))
 }
\seealso{\code{\link{distributionSet}}}
\keyword{models}
\keyword{distribution}
