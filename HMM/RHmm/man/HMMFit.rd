\name{HMMFit}
\alias{HMMFit}
\alias{HMMFitClass}
\alias{summary.HMMFitClass}
\alias{print.summary.HMMFitClass}
\title{Fit an Hidden Markov Model}
\description{This function returns an HMMFitClass object which contains the results
of the Baum-Welch algorithm for the user's data}
\usage{
   HMMFit(obs, dis="NORMAL", nStates=2, asymptCov=FALSE, ...)
}
\details{
Typical usages are:
\itemize{
    \item HMMFit(obs, dis="NORMAL", nStates=, ..., asymptCov=FALSE)
    \item HMMFit(obs, dis="DISCRETE", nStates=, levels=NULL, ..., asymptCov=FALSE)
    \item HMMFit(obs, dis="MIXTURE", nStates=, nMixt=, ..., asymptCov=FALSE)
}
The differents arguments are:
    \describe{
    \item{obs}{A vector, a matrix, a data frame, a list of vectors or a list of matrices of observations. See section \bold{obs parameter}.}
    \item{dis}{Distribution name. In 'NORMAL', 'DISCRETE' or 'MIXTURE'. Default 'NORMAL'.}
    \item{nStates}{Number of hidden states. Default 2.}
    \item{nMixt}{Number of mixtures of normal distributions if dis ='MIXTURE'}
    \item{levels}{A character vector of all different levels of 'obs'. By Default (levels=NULL), this vector is computed from 'obs'.}
    \item{asymptCov}{A boolean. asymptCov=TRUE if the asymptotic covariance matrix is computed. Default FALSE.}
    \item{...}{optional parameter:\describe{
        \item{control: }{A list of control parameters for the Baum-Welch algorithm. See \bold{control parameter}}
        }}
    }
}
\arguments{
    \item{obs}{A vector, a matrix, a data frame, a list of vectors or a list of matrices of observations. See section \bold{obs parameter}.}
    \item{dis}{Distribution name. In 'NORMAL', 'DISCRETE' or 'MIXTURE'. Default 'NORMAL'.}
    \item{nStates}{Number of hidden states. Default 2.}
     \item{asymptCov}{A boolean. asymptCov=TRUE if the asymptotic covariance matrix is computed. Default FALSE.}
    \item{...}{Other parameters. See \bold{Details}.}
        
}
\value{ a \sQuote{HMMFitClass} object with the following elements:
    \describe{
        \item{HMM}{A HMMClass object with the fitted values of the model. See \bold{HMMSet.}}
        \item{LLH}{Log-likelihood}
        \item{BIC}{BIC criterium}
        \item{nIter}{Number of iterations of the Baum-Welch algorithm}
        \item{relVariation}{last relative variation of the LLH function}
        \item{asymptCov}{Asymptotic covariance matrix of independant parameters. NULL if not computed.}
        \item{obs}{the observations.}
        \item{call}{The call object of the function call}
        }
    }
\section{obs parameter}{
        If you fit the model with only one sample, obs is
        either a vector (for univariate distributions) or a matrix (for multivariate distributions) or a data frame.
        In the two last cases, the number of columns of obs defines the dimension of observations.\cr\cr

        If you fit the model with more than one sample, obs is a list of samples. Each element of obs is then a vector
        (for univariate distributions) or a matrix (for multivariate distributions). The samples do not need to have the same length.\cr\cr

        For discrete distributions, obs can be a vector (or a list of vectors) of any type of R factor objects.
    }


\section{control parameter}{
    \describe{
        \item{init}{Kind of initialisation ='KMEANS' (for univariate or multivariate normal distributions), 'RANDOM' or 'USER'. Default 'RANDOM', see \bold{Random Initialization}}
        \item{iter}{Maximum number of iterations for the Baum-Welch algorithm. Default 500}
        \item{tol}{Tolerance of the relative log-likehood augmentation. Default 1e-6}
        \item{verbose}{=0, no details, =1 iterations are displayed. Default 0}
        \item{nInit}{Number of random initialisations. Default 5}
        \item{nIterInit}{Number of maximum iterations of the Baum-Welch algorithm in the random initialisation phase. Default 5}
        \item{initPoint}{An HMMClass object used to initialize the parameters of the Baum-Welch algorithm. Default NULL.\cr
            if initPoint != NULL, init is set to "USER"}
        }
}

\section{random initialization}{
    'initProb' and 'transMat' parameters are uniformly drawn.\cr\cr
    For univariate normal distributions, empirical mean \eqn{m}{m} and variance \eqn{\sigma^2}{s^2}
    of all the samples are computed.
    Then for every states,
    an initial value of the 'mean' parameter is uniformly drawn between \eqn{m - 3\,\sigma}{m - 3s} and \eqn{m - 3\,\sigma}{m + 3s}
    and an initial value of the 'var'
    parameter is uniformly drawn between \eqn{\frac{1}{2}\,\sigma^2}{0.5 s^2} and \eqn{3\,\sigma^2}{3 s^2}.
    \cr
    For multivariate normal distributions, the same procedure is applied for each component of the mean vectors.
    The initial covariance matrix is diagonal, and each initial variance is computed as for univariate models.
    \cr
    For mixtures of univariate normal distributions, initial values for 'mean' and 'var' parameters are computed
    the same way than for normal distributions. The initial value of 'proportion' parameter is uniformly drawn.
    \cr
     For mixtures of multivariate normal distributions, the same procedure is applied for each component of the mean vectors,
     all the covariance matrices are diagonal and each initial variance is computed as for univariate models. The initial value
     of 'proportion' parameter is also uniformly drawn.
    \cr
    For discrete distributions, the initial values of 'proba' parameters are uniformly drawn.
    \cr
    Of course, the initial values of the parameters 'initProba', 'proba', 'proportion' and 'transMat' are standardized to
    ensure that they can represent probabilities vectors or transition matrices.
}

\section{asymptotic covariance matrix}{
    The asymptotic covariance matrix of estimates is computed using the Lystig and Hugues's algorithm.See \bold{asymptoticCov}.
    \cr
    The summary and print.summary methods display the results.
}

\references{
    Bilmes Jeff A. (1997) \emph{ A Gentle Tutorial of the EM Algorithm and its Application to Parameter
    Estimation for Gaussian Mixture and Hidden Markov Models} \url{http://ssli.ee.washington.edu/people/bilmes/mypapers/em.ps.gz}

 }

\examples{
    # Fit a 3 states 1D-gaussian model
    data(n1d_3s)
    HMMFit(obs_n1d_3s, nStates=3)

    # Fit a 3 states gaussian HMM for obs_n1d_3s
    # with iterations printing and kmeans initialization
    Res_n1d_3s <- HMMFit(obs=obs_n1d_3s, nStates=3,
        control=list(verbose=1, init="KMEANS"),
        asymptCov=TRUE)
    summary(Res_n1d_3s)

    # Fit a 2 states 3D-gaussian model
    data(n3d_2s)
    summary(HMMFit(obs_n3d_2s, asymptCov=TRUE))

    # Fit a 2 states mixture of 3 normal distributions HMM
    # for data_mixture
    data(data_mixture)
    ResMixture <- HMMFit(data_mixture, nStates=2, nMixt=3,
        dis="MIXTURE")

    # Fit a 3 states discrete HMM for weather data
    data(weather)
    ResWeather <- HMMFit(weather, dis='DISCRETE', nStates=3)
 }
\keyword{htest}
