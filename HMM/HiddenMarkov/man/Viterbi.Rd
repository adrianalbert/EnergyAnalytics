\name{Viterbi}
\alias{Viterbi}
\alias{Viterbi.dthmm}
\alias{Viterbi.mmglm0}
\alias{Viterbi.mmglm1}
\alias{Viterbi.mmglmlong1}

\title{Viterbi Algorithm for Hidden Markov Model}

\description{
Provides methods for the generic function \code{\link{Viterbi}}. This predicts the most likely sequence of Markov states given the observed dataset. There is currently no method for objects of class \code{"\link{mmpp}"}.
}

\usage{
\method{Viterbi}{dthmm}(object, ...)
\method{Viterbi}{mmglm0}(object, ...)
\method{Viterbi}{mmglm1}(object, ...)
\method{Viterbi}{mmglmlong1}(object, ...)
}

\arguments{
  \item{object}{an object with class \code{"\link{dthmm}"}, \code{"\link{mmglm0}"}, \code{"\link{mmglm1}"} or \code{"\link{mmglmlong1}"}.}
  \item{...}{other arguments.}
}

\details{
The purpose of the Viterbi algorithm is to \emph{globally decode} the underlying hidden Markov state at each time point. It does this by determining the sequence of states \eqn{(k_1^*, \cdots, k_n^*)}{(k_1^*, ..., k_n^*)} which maximises the joint distribution of the hidden states given the entire observed process, i.e.
\deqn{
(k_1^*, \cdots, k_n^*) = _{\stackrel{\mbox{argmax}}{{k_1, \cdots, k_n \in \{1, 2, \cdots, m\}}}} \Pr\{ C_1=k_1, \cdots, C_n=k_n \,|\, X^{(n)}=x^{(n)} \}\,.
}{
(k_1^*, ..., k_n^*) = argmax Pr{ C_1=k_1, ..., C_n=k_n | X^{(n)}=x^{(n)} } .
}

The algorithm has been taken from Zucchini (2005), however, we calculate sums of the logarithms of probabilities rather than products of probabilities. This lessens the chance of numerical underflow. Given that the logarithmic function is monotonically increasing, the \emph{argmax} will still be the same. Note that \emph{argmax} can be evaluated with the \R function \code{\link[base:which.min]{which.max}}.

Determining the \emph{a posteriori} most probable state at time \eqn{i} is referred to as \emph{local decoding}, i.e.
\deqn{
k_i^* = _{\stackrel{\mbox{argmax}}{k \in \{1, 2, \cdots, m\}}} \Pr\{ C_i=k \,|\, X^{(n)}=x^{(n)} \}\,.
}{
k_i^* = argmax Pr{ C_i=k | X^{(n)}=x^{(n)} } .
}
Note that the above probabilities are calculated by the function \code{\link{Estep}}, and are contained in \code{u[i,j]} (output from \code{Estep}), i.e. \eqn{k_i^*} is simply \code{which.max(u[i,])}.

The code for the methods \code{"\link{dthmm}"}, \code{"\link{mmglm0}"}, \code{"\link{mmglm1}"} or \code{"\link{mmglmlong1}"} can be viewed by typing \code{Viterbi.dthmm}, \code{Viterbi.mmglm0}, \code{Viterbi.mmglm1} or \code{Viterbi.mmglmlong1}, respectively, on the \R command line.
}

\value{
A vector of length \eqn{n} containing integers (\eqn{1, \cdots, m}{1, ..., m}) representing the hidden Markov states for each node of the chain.
}

\references{
Cited references are listed on the \link{HiddenMarkov} manual page.
}

\examples{
Pi <- matrix(c(1/2, 1/2,   0,   0,   0,
               1/3, 1/3, 1/3,   0,   0,
                 0, 1/3, 1/3, 1/3,   0,
                 0,   0, 1/3, 1/3, 1/3,
                 0,   0,   0, 1/2, 1/2),
             byrow=TRUE, nrow=5)
delta <- c(0, 1, 0, 0, 0)
lambda <- c(1, 4, 2, 5, 3)
m <- nrow(Pi)

x <- dthmm(NULL, Pi, delta, "pois", list(lambda=lambda), discrete=TRUE)
x <- simulate(x, nsim=2000)

#------  Global Decoding  ------

states <- Viterbi(x)
states <- factor(states, levels=1:m)

#  Compare predicted states with true states
#  p[j,k] = Pr{Viterbi predicts state k | true state is j}
p <- matrix(NA, nrow=m, ncol=m)
for (j in 1:m){
    a <- (x$y==j)
    p[j,] <- table(states[a])/sum(a)
}
print(p)

#------  Local Decoding  ------

#   locally decode at i=100

print(which.max(Estep(x$x, Pi, delta, "pois", list(lambda=lambda))$u[100,]))

#---------------------------------------------------
#   simulate a beta HMM

Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)
delta <- c(0, 1)

y <- seq(0.01, 0.99, 0.01)
plot(y, dbeta(y, 2, 6), type="l", ylab="Density", col="blue")
points(y, dbeta(y, 6, 2), type="l", col="red")

n <- 100
x <- dthmm(NULL, Pi, delta, "beta",
           list(shape1=c(2, 6), shape2=c(6, 2)))
x <- simulate(x, nsim=n)

#   colour denotes actual hidden Markov state
plot(1:n, x$x, type="l", xlab="Time", ylab="Observed Process")
points((1:n)[x$y==1], x$x[x$y==1], col="blue", pch=15)
points((1:n)[x$y==2], x$x[x$y==2], col="red", pch=15)

states <- Viterbi(x)
#   mark the wrongly predicted states
wrong <- (states != x$y)
points((1:n)[wrong], x$x[wrong], pch=1, cex=2.5, lwd=2)
}

\keyword{methods}
\keyword{optimize}
