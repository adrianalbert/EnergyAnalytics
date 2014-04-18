\name{HMMGraphicDiag}
\alias{HMMGraphicDiag}
\title{Graphic diagnostic of the HMM estimation}
\description{This function plots the kernel density of the observations and the normal (mixture of normal, discrete) density
    with estimated parameters for each hidden states.}
\usage{
HMMGraphicDiag(vit, HMM, obs, color="green")
}
\arguments{
    \item{vit}{A ViterbiClass object which gives the hidden states}
    \item{HMM}{A HMMClass or a HMMFitClass object which describes the model}
    \item{obs}{The vector, list of vectors of observations}
    \item{color}{Color for the kernel density plot}
}

\value{None.}
\note{
 HMMGraphicDiag is not implemented for multivariate distributions.\cr

 The kernel densities of observations for each hidden states of the model are plotting using:\cr
    plot(density(obs[vit$states=i])) and i in 1..HMM$nStates (or HMM$HMM$nStates)

  }
\examples{
data(n1d_3s)
obs <- obs_n1d_3s
#Fits an 3 states gaussian model
ResFit <- HMMFit(obs, nStates=3)
VitPath <- viterbi(ResFit, obs)
# Graphic diagnostic
HMMGraphicDiag(VitPath, ResFit, obs)
}
\seealso{HMMFit, viterbi}
\keyword{htest}
