\name{predict.protoclass}
\alias{predict.protoclass}
\title{Make predictions based on prototypes}
\usage{
  \method{predict}{protoclass} (object, newx, z, ...)
}
\arguments{
  \item{object}{output of protoclass function}

  \item{newx}{features from some new data that we want to
  make predictions on.}

  \item{z}{potential prototype positions z (often taken to
  be x, the training data)}

  \item{...}{ignore this argument.}
}
\description{
  Performs nearest-prototype classification on newx.
}
\seealso{
  \code{\link{protoclass}}
  \code{\link{predictwithd.protoclass}}
}

