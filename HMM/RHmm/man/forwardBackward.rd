\name{forwardBackward}
\alias{forwardbackward}
\alias{forwardBackward}
\title{forward-backward function}
\description{The forward-backward function is used to compute quantities used in the Baum-Welch algorithm.}
\usage{
forwardBackward(HMM, obs, logData=TRUE)
}
\arguments{
    \item{HMM}{a HMMClass or a HMMFitClass object}
    \item{obs}{a vector (matrix) of observations, or a list of vectors (or matrices) if there are more than one samples}
    \item{logData}{a boolean. If true, the function computes the logaritm of the Alpha, Beta and Rho quantities instead of the quantities themselves.}
    }

\value{ If obs is one sample, a list of following elements, if obs is a list of samples, a
    list of list of following elements. See \bold{note} for mathematical definitions.
    \item{Alpha}{The matrix of (log) 'forward' probabilities (density) (size: number of obs. times number of hidden states)}
    \item{Beta}{The matrix of (log) 'backward' probabilities (density) (size: number of obs. times number of hidden states)}
    \item{Delta}{The matrix of 'conditional forward' probabilities (density) (size: number of obs. times number of hidden states)}
    \item{Gamma}{The matrix of probabilities (density) of being at time t in state i (size: number of obs. times number of hidden states)}
    \item{Xsi}{The matrix of probabilities (density) of being in state i at time t and being in state j at time t + 1 (size: number of obs. times number of hidden states)}
    \item{Rho}{The vector of (log) probabilities (density) of seeing the partial sequence obs[1] \ldots obs[t] (size number of obs.)}
    \item{LLH}{Log-likelihood}
     }

\note{
 Let \eqn{o=(o(1),\,\ldots,\,o(T))}{obs=(obs[1], \ldots obs[T])} be the
 vector of observations, and \eqn{O=(O(t), t=1,\,\ldots,\,T)}{O=(O[t],
 t=1, \ldots, T)}, the corresponding random variables. Let \eqn{Q=(Q(t), t=1,\,\ldots,\,T)}{(Q[t], t=1, \ldots, T)}
 be the hidden Markov chain whose values are in \eqn{\left\{1,\,\ldots,\,nStates\right\}}{{1, \ldots, nStates}}
 We have the
 following definitions:\cr

 \eqn{\alpha_i(t) =
 P(O_1=o(1),\,\ldots,\,O(t)=o(t),\,Q(t)=i\,|\,HMM)}{Alpha[i][t] =
 P(O[1]=obs[1], \ldots, O[t]=obs[t], Q[t]=i | HMM)} which is
 the probability of seeing the partial sequence
 \eqn{o(1),\,\ldots,\,o(t)}{obs[1], \ldots, obs[t]} and ending up
 in state i at time t.\cr

 \eqn{\beta_i(t) = P(O_{t+1}=o(t+1),\,\ldots,\,O(T)=o(T),\,Q(t)=i
| HMM)}{Beta[i][t] =
 P(O[t+1]=obs[t+1], \ldots, O[T]=obs[T], Q[t]=i | HMM)} which
 is the probability of the ending partial sequence \eqn{o(t+1),\,\ldots,\,o(T)}{obs[t+1], \ldots, obs[T]}
 given that we started at state i at time t.\cr

 \eqn{\delta_i(t) =
 P(Q(t)=i\,|\,O_1=o(1),\,\ldots,\,O(t)=o(t),\,HMM)}{Delta[i][t] =
 P(Q[t]=i | O[1]=obs[1], \ldots, O[t]=obs[t], HMM)} which is
 the probability of beeing in state i at time t knowing the  partial sequence
 \eqn{o(1),\,\ldots,\,o(t)}{obs[1], \ldots, obs[t]}.\cr

 \eqn{\Gamma_i(t) = P(Q(t) = i\,|\,O=o,\,HMM)}{Gamma[i][t] = P(Q[t]=i | O=obs, HMM)} which is the probability of being in state i
 at time t for the state sequence \eqn{O=o}{O=obs}. \cr
 \eqn{\xi_{ij}(t)=P(Q(t)=i,\,Q(t+1)=j\,|\,O=o,\,HMM)}{Xsi[i][j][t]=P(Q[t]=i, Q[t+1]=j | O=obs, HMM)} which is the probability of being
 in state i at time t and being in state j at time t + 1.\cr

 \eqn{\rho(t) = P(O_1=o(1),\,\ldots,\,O_t=o(t)\,|\, HMM)}{Rho[t] = P(O[1]=obs[1], \ldots, O[t]=obs(t) | HMM)} witch is probabilities of seeing
 the partial sequence \eqn{o(1),\,\ldots,\,o(t)}{obs[1] \ldots obs[t]}.\cr

 \eqn{LLH=\ln\rho[T]}{LLH=ln(Rho[T])}
 \cr
  When the sequences of observations become larger, the probabilistic values in this algorithm get increasingly small and after enough iterations become almost zero.
  For that reason, the Alphas, Betas and Rhos are scaled during the iterations of the algorithm to avoid underflow problems. The logarithm of these probabilistic values are compute
  from the logarithm of the scaled quantities and should produce a more precise result.
  \cr
  If the conditional distribution is continuous, the probabilistic values are replaced by density values and can be greater than one.

}
\references{
    Jeff A. Bilmes (1997) \emph{ A Gentle Tutorial of the EM Algorithm and its Application to Parameter
    Estimation for Gaussian Mixture and Hidden Markov Models} \url{http://ssli.ee.washington.edu/people/bilmes/mypapers/em.ps.gz}
}
\examples{
    data(n1d_3s)
    #Fits an 2 states gaussian model for geyser duration
    Res_n1d_3s <- HMMFit(obs_n1d_3s, nStates=3)
    #Forward-backward procedure with log Alpha, Beta, Rho
    fbLog <- forwardBackward(Res_n1d_3s, obs_n1d_3s)
    #Forward-backward procedure with Alpha Beta and Rho
    fb <- forwardBackward(Res_n1d_3s, obs_n1d_3s, FALSE)
  }

\keyword{htest}
