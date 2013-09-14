\name{predictwithd.protoclass}
\alias{predictwithd.protoclass}
\title{Make predictions based on prototypes}
\usage{
  predictwithd.protoclass(object, d.test.z)
}
\arguments{
  \item{object}{output of protoclass function}

  \item{d.test.z}{ntest-by-m matrix of distances between
  test features and z}
}
\description{
  Performs nearest-prototype classification like
  predict.protoclass, except only requires the ntest-by-m
  matrix of distances between test features and z rather
  than requiring these matrices themselves. (Using this
  saves repeated computation of the pairwise distances...)
}
\seealso{
  \code{\link{predict.protoclass}}
}

