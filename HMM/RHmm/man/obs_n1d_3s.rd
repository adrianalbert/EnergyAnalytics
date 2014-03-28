\name{obs_n1d_3s}
\docType{data}
\alias{n1d_3s}
\alias{obs_n1d_3s}
\title{A 3 states HMM with univariate gaussian distribution data set}
\description{This data set contains 10 samples with 100 observations of this HMM.}
\details{
    File 'n1d_3s.RData' contains the list of vectors 'obs_n1d_3s'.
}
\note{Parameters for the simulation
\cr
Model:\cr
------\cr
3 states HMM with univariate gaussian distribution\cr
\cr
Initial probabilities:\cr
  \tabular{rrr}{
  Pi1 \tab Pi2 \tab Pi3\cr
  0.2 \tab 0.4 \tab 0.4\cr
  }

Transition matrix:\cr
 \tabular{crrr}{
        \tab State 1 \tab State 2 \tab State 3\cr
State 1 \tab 0.5 \tab 0.1 \tab 0.4\cr
State 2 \tab 0.2 \tab 0.7 \tab 0.1\cr
State 3 \tab 0.3 \tab 0.3 \tab 0.4\cr
}

Conditionnal distribution parameters:\cr
\cr
Distribution parameters:\cr
 \tabular{crr}{
     \tab mean \tab var\cr
State 1 \tab 10 \tab 4\cr
State 2 \tab -5 \tab 2\cr
State 3 \tab -1 \tab 1\cr
	}
}
\usage{data(n1d_3s)}
\format{RData file}

\keyword{datasets}
