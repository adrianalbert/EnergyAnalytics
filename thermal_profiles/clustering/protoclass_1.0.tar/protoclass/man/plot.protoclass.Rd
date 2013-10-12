\name{plot.protoclass}
\alias{plot.protoclass}
\title{Plots prototypes with data}
\usage{
  \method{plot}{protoclass} (x, xtr, y, z = NULL,
    circleFLAG = TRUE, boundFLAG = FALSE, res = NULL,
    main = NULL, ...)
}
\arguments{
  \item{x}{output of protoclass function}

  \item{xtr}{n by p matrix of training features.}

  \item{y}{n-vector of labels of the training data.}

  \item{z}{set of potential prototypes (only required if
  protoclass was not run with z=x).}

  \item{circleFLAG}{indicates whether to draw circles
  around each prototype.}

  \item{boundFLAG}{indicates whether to draw classification
  boundary.  This can be computationally expensive.}

  \item{res}{resolution of the classification boundary.
  Computation increases with res.  Try 100.}

  \item{main}{optional main title for plot.}

  \item{...}{other arguments to pass to plot.}
}
\description{
  This function is only for two-dimensional data.
}
\seealso{
  \code{\link{protoclass}}
}

