\name{data_mixture}
\docType{data}
\alias{data_mixture}
\title{Simulated univariate mixture of 3 gausssian distributions}
\description{This data set contains 10 samples with 100 observations of this HMM}
\note{Parameters for the simulation
Model:\cr
------\cr
2 states HMM with gaussian mixture distribution\cr
\cr
Initial probabilities:\cr
  \tabular{rr}{
 Pi1 \tab Pi2\cr
  0.4 \tab 0.6\cr
}
Transition matrix:\cr
   \tabular{crr}{
      \tab State 1 \tab State 2\cr
State 1 \tab 0.8 \tab 0.2\cr
State 2 \tab 0.4 \tab 0.6\cr
}
Conditionnal distribution parameters:\cr

State 1\cr
\tabular{crrr}{ 
\tab mean \tab var \tab prop\cr
mixt. 1 \tab 1 \tab 2 \tab 0.1\cr
mixt. 2 \tab 2 \tab 3 \tab 0.2\cr
mixt. 3 \tab 3 \tab 4 \tab 0.7\cr
}
State 2\cr
\tabular{crrr}{ 
\tab mean \tab var \tab prop\cr
mixt. 1 \tab 1 -1 \tab 1 4 \tab 1 0.4\cr
mixt. 2 \tab 1 -2 \tab 1 2 \tab 1 0.3\cr
mixt. 3 \tab 1 -3 \tab 1 1 \tab 1 0.3\cr
}
}

\usage{data(data_mixture)}
\format{R script}

\keyword{datasets}
