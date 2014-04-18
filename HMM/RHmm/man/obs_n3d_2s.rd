\name{obs_n3d_2s}
\docType{data}
\alias{n3d_2s}
\alias{obs_n3d_2s}
\title{A 2 states HMM with 3D gaussian distribution data set}
\description{This data set contains 1 sample with 1000 observations of this HMM.}
\details{
    File 'n3d_3s.RData' contains the matrix 'obs_n3d_2s'.
}

\note{Parameters for the simulation
\cr
Model:\cr
------\cr
2 states HMM with 3-d gaussian distribution\cr
\cr
Initial probabilities:\cr
 \tabular{rr}{
   Pi1 \tab Pi2\cr
  0.2 \tab 0.8\cr
  }

Transition matrix:\cr
    \tabular{crr}{
    \tab State 1 \tab State 2\cr
    State 1 \tab 0.1 \tab 0.9\cr
    State 2 \tab 0.8\tab 0.2\cr
}

Conditionnal distribution parameters:\cr
\cr
Distribution parameters:\cr
  State 1\cr
    \tabular{rrrr}{
     mean \tab \tab \tab cov matrix\cr
        2 \tab 1.0 \tab 1.6 \tab -0.6\cr
        2 \tab 1.6 \tab 4.0 \tab -3.6\cr
        2 \tab -0.6 \tab -3.6 \tab 9.0\cr
        }

  State 2\cr
   \tabular{rrrr}{
    mean \tab cov matrix \tab \tab \cr
    -1 \tab 4.0 \tab 1.6 \tab 1.6\cr
    -2 \tab 1.6 \tab 1.0 \tab 1.0\cr
    -3 \tab 1.6 \tab 1.0 \tab 4.0\cr
    }
}

\usage{data(n3d_2s)}
\format{RData file}

\keyword{datasets}
