\name{HiddenMarkov-internal}
\alias{residuals.mmpp}
\alias{Viterbi.mmpp}
\alias{as.dthmm}
\alias{as.mmglm0}
\alias{makedensity}
\alias{makedensity1}
\alias{makedistn}
\alias{getj}
\alias{dglm}
\alias{pglm}
\alias{dmmglm}
\alias{pmmglm}
\alias{mmglm}

\title{Internally Used Functions}

\description{
This page lists internally used functions. They should not be required by most users.
}

\usage{
as.dthmm(object)
as.mmglm0(object)

makedensity(distn)
makedensity1(distn)
makedistn(distn)

getj(x, j)

dglm(x, x1, beta0, beta1, sigma, family, link, size = NA,
     log = FALSE)
pglm(q, x1, beta0, beta1, sigma, family, link, size = NA,
     log.p = FALSE)
dmmglm(x, beta, sigma, glmfamily, Xdesign, size=NA, log=FALSE)
pmmglm(x, beta, sigma, glmfamily, Xdesign, size=NA, log.p=FALSE)

mmglm(...)
}

\details{
The function \code{makedensity} is used to reparameterise various \R probability (or density) functions (e.g. \code{\link{dnorm}} and \code{\link{dpois}}) into a format with a standard argument list. Similarly, \code{makedistn} reparameterises a distribution function.

The function \code{getj} is used to extract the \eqn{j}th element from each vector component in a list object.

The function \code{as.dthmm} coerces an object with class \code{"\link{mmglm0}"} to an object with class \code{"\link{dthmm}"}. Similarly, the function \code{as.mmglm0} coerces an object with class \code{"\link{dthmm}"} to an object with class \code{"\link{mmglm0}"} (if possible).

The functions \code{dglm} and \code{pglm} calculate the density and probability, respectively, for an observation given a generalised linear model.

The functions \code{dmmglm} and \code{pmmglm} calculate the density and probability, respectively, for an observation given a generalised linear model.
}

\keyword{internal}
