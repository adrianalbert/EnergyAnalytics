\name{HMMSim}
\alias{HMMSim}
\title{Simulation of an Hidden Markov Model}
\description{Simulation of an HMM for different classes of observations distributions.}
\usage{
HMMSim(nSim, HMM, lastState=NULL)
}
\arguments{
    \item{nSim}{Number of simulations.}
    \item{HMM}{An HMMClass object. See \bold{HMMSet}}
    \item{lastState}{If not NULL, the state of the previous observation (usefull to complete a serie)}
    }
\value{ a list with
    \item{obs}{simulated observations (a vector for univariate distributions, a matrix for multivariate distributions)}
    \item{states}{simulated hidden states}
    }
\examples{
    # simulate a 3 hidden states model with univariate normal distributions
    n_1d_3s <- distributionSet("NORMAL", mean=c(1, -2, 5), var=c(1, 2, 4))
    initProb3 <- rep(1,3)/3
    transMat3 <- rbind(c(0.5, 0.4, 0.1), c(0.3, 0.4, 0.3), c(0.2, 0.1, 0.7))
    hmm_1d_3s <- HMMSet(initProb3, transMat3, n_1d_3s)
    simul <- HMMSim(1000, hmm_1d_3s)
    #Complete the serie.
    simul <- c(simul, HMMSim(1000, hmm_1d_3s, simul$states[1000]))
    }
 \seealso{\code{\link{HMMSet}}}
\keyword{datagen}
\keyword{distribution}
