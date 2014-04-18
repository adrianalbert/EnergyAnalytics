\name{bwcontrol}
\alias{bwcontrol}

\title{Control Parameters for Baum Welch Algorithm}

\description{
Creates a list of parameters that control the operation of \code{\link{BaumWelch}}.
}

\usage{
bwcontrol(maxiter = 500, tol = 1e-05, prt = TRUE, posdiff = TRUE,
          converge = expression(diff < tol))
}

\arguments{
  \item{maxiter}{is the maximum number of iterations, default is 500.}
  \item{tol}{is the convergence criterion, default is 0.00001.}
  \item{prt}{is logical, and determines whether information is printed at each iteration; default is \code{TRUE}.}
  \item{posdiff}{is logical, and determines whether the iterative process stops if a negative log-likelihood difference occurs, default is \code{TRUE}.}
  \item{converge}{is an expression giving the convergence criterion. The default is the difference between successive values of the log-likelihood.}
}

\examples{
#    Increase the maximum number of iterations to 1000.
#    All other components will retain their default values.
a <- bwcontrol(maxiter=1000)
print(a)
}

\keyword{optimize}
