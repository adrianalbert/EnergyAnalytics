\name{BaumWelch}
\alias{BaumWelch}
\alias{BaumWelch.dthmm}
\alias{BaumWelch.mmglm0}
\alias{BaumWelch.mmglm1}
\alias{BaumWelch.mmglmlong1}
\alias{BaumWelch.mmpp}

\title{Estimation Using Baum-Welch Algorithm}

\description{
Estimates the parameters of a hidden Markov model. The Baum-Welch algorithm (Baum et al, 1970) referred to in the HMM literature is a version of the EM algorithm (Dempster et al, 1977). See Hartley (1958) for an earlier application of the EM methodology, though not referred to as such.}

\usage{
BaumWelch(object, control, ...)
\method{BaumWelch}{dthmm}(object, control = bwcontrol(), ...)
\method{BaumWelch}{mmglm0}(object, control = bwcontrol(), ...)
\method{BaumWelch}{mmglm1}(object, control = bwcontrol(), ...)
\method{BaumWelch}{mmglmlong1}(object, control = bwcontrol(), SNOWcluster=NULL,
          tmpfile=NULL, ...)
\method{BaumWelch}{mmpp}(object, control = bwcontrol(), ...)
}

\arguments{
  \item{object}{an object of class \code{"\link{dthmm}"}, \code{"\link{mmglm0}"}, \code{"\link{mmglm1}"}, \code{"\link{mmglmlong1}"}, or \code{"\link{mmpp}"}.}
  \item{control}{a list of control settings for the iterative process. These can be changed by using the function \code{\link{bwcontrol}}.}
  \item{SNOWcluster}{see section below called \dQuote{Parallel Processing}.}
  \item{tmpfile}{name of a file (.Rda) into which estimates are written at each 10th iteration. The model object is called \code{object}. If \code{NULL} (default), no file is created.}
  \item{\dots}{other arguments.}
}

\details{
The initial parameter values used by the EM algorithm are those that are contained within the input \code{object}.

The code for the methods \code{"\link{dthmm}"}, \code{"\link{mmglm0}"}, \code{"\link{mmglm1}"},\code{"\link{mmglmlong1}"},  and \code{"\link{mmpp}"} can be viewed by typing \code{BaumWelch.dthmm}, \code{BaumWelch.mmglm0}, \code{BaumWelch.mmglm1}, \code{BaumWelch.mmglmlong1}  or \code{BaumWelch.mmpp}, respectively, on the \R command line.
}

\value{
The output object (a \code{\link[base]{list}}) with have the same class as the input, and will have the same components. The parameter values will be replaced by those estimated by this function. The object will also contain additional components. 

An object of class \code{"\link{dthmm}"} will also contain
  \item{u}{an \eqn{n \times m}{n*m} matrix containing estimates of the conditional expectations. See \dQuote{Details} in \code{\link{Estep}}.}
  \item{v}{an \eqn{n \times m \times m}{n*m*m} array containing estimates of the conditional expectations. See \dQuote{Details} in \code{\link{Estep}}.}
  \item{LL}{value of log-likelihood at the end.}
  \item{iter}{number of iterations performed.}
  \item{diff}{difference between final and previous log-likelihood.}
}

\section{Parallel Processing}{
In longitudinal models, the forward and backward equations need to be calculated for each individual subject. These can be done independently, the results being concatenated to be used in the E-step. If the argument \code{SNOWcluster} is set, subjects are divided equally between each node in the cluster for the calculation of the forward and backward equations. This division is very basic, and assumes that all nodes run at a roughly comparable speed.

If the communication between nodes is slow and the dataset is small, then the time taken to allocate the work to the various nodes may in fact take more time than simply using one processor to perform all of the calculations.

The required steps in initiating parallel processing are as follows.
\preformatted{
#   load the "snow" package
library(snow)

#   define the SNOW cluster object, e.g. a SOCK cluster
#   where each node has the same R installation.
cl <- makeSOCKcluster(c("localhost", "horoeka.localdomain", 
                        "horoeka.localdomain", "localhost"))

#   A more general setup: Totara is Fedora, Rimu is Debian:
#   Use 2 processors on Totara, 1 on Rimu:
totara  <- list(host="localhost",
                rscript="/usr/lib/R/bin/Rscript",
                snowlib="/usr/lib/R/library")
rimu    <- list(host="rimu.localdomain",
                rscript="/usr/lib/R/bin/Rscript",
                snowlib="/usr/local/lib/R/site-library")
cl <- makeCluster(list(totara, totara, rimu), type="SOCK")

#   then define the required model object
#   say the model object is called x
BaumWelch(x, SNOWcluster=cl)

#   stop the R jobs on the slave machines
stopCluster(cl)
}
Note that the communication method does not need to be \code{SOCKS}; see the \pkg{snow} package documentation, topic \code{\link[snow]{snow-startstop}}, for other options. Further, if some nodes are on other machines, the firewalls may need to be tweaked. The master machine initiates the \R jobs on the slave machines by communicating through port 22 (use of security keys are needed rather than passwords), and subsequent communications through port 10187. Again, these details can be tweaked in the options settings within the \pkg{snow} package.
}

\seealso{
\code{\link{logLik}}, \code{\link{residuals}}, \code{\link{simulate}}, \code{\link{summary}}, \code{\link{neglogLik}}}

\references{
Cited references are listed on the \link{HiddenMarkov} manual page.
}

\keyword{methods}
\keyword{optimize}
