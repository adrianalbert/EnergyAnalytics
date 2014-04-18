\name{mmglm-2nd-level-functions}
\alias{mmglm-2nd-level-functions}
\alias{Estep.mmglm1}
\alias{Mstep.mmglm1}

\title{Markov Modulated Generalised Linear Model - 2nd Level Functions}

\description{
These functions will generally not be called directly by the user.
}

\usage{
Estep.mmglm1(object, fortran=TRUE)
Mstep.mmglm1(object, u)
}

\arguments{
  \item{object}{an object with class \code{"mmglm"} or \code{"mmglmlong"}}
  \item{u}{a matrix of weights by Markov state and observation used in fitting the generalised linear model.}
  \item{fortran}{logical, if \code{TRUE} (default) use the Fortran code in the forward-backward equations, else use the \R code.}
}

\keyword{misc}

