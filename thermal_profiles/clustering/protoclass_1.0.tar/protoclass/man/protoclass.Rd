\docType{package}
\name{protoclass}
\alias{protoclass}
\alias{protoclass-package}
\title{Greedy algorithm for prototype selection}
\usage{
  protoclass(x, y, z, dxz, eps, lambda = 1/n)
}
\arguments{
  \item{x}{n by p matrix of training features (optional,
  see \code{dxz}).}

  \item{y}{n-vector of labels of the training data.}

  \item{z}{set of potential prototypes (optional, see
  \code{dxz}).}

  \item{dxz}{instead of \code{x} and \code{z}, you can give
  \code{dxz}, the matrix of pairwise dissimilarities
  between \code{x} and \code{z}, with ij-th element giving
  the dissimilarity between training point x_i and
  prototype-candidate z_j.}

  \item{eps}{size of covering balls.}

  \item{lambda}{cost of adding a prototype.}
}
\value{
  An object of class "protoclass," which has the following
  elements: \itemize{ \item{alpha: }{Matrix of dimensions
  nrow(z)-by-nclass.  alpha[j,k] indicates whether jth
  potential prototype has been chosen as a prototype for
  class k.} \item{classes: }{Names of classes}
  \item{proto.order: }{The sequence of prototypes that were
  selected.} \item{ncovered: }{nproto-by-nclass matrix with
  ncovered[j,k] giving the number of class k training
  points covered by the jth prototype's ball.}
  \item{coverlist: }{n-by-nclass matrix with row i giving
  number of each type of prototype covering point i.}
  \item{uncovered: }{Indicates whether a training point is
  not covered by a prototype of its own class.}
  \item{wrongcover: }{Number of prototypes from other
  classes covering each training point.} \item{nproto:
  }{nclass-vector giving the number of prototypes in each
  class.} }
}
\description{
  Selects prototypes for each class in a greedy manner as
  described in 'Bien and Tibshirani (2011) Prototype
  Selection for Interpretable Classification. Annals of
  Applied Statistics. 5(4). 2403-2424.'

  protoclass
}
\details{
  It's more efficient to compute dxz just once on your own
  rather than have protoclass repeatedly compute the
  pairwise distances on each call.
}
\examples{
# generate some data:
set.seed(1)
n <- 200
p <- 2
x <- matrix(rnorm(n * p), n, p)
y <- rep(c("A","B"), each=n/2)
x[y=="A", ] <- x[y=="A", ] + 3
itr <- sample(n, n/2)
xtr <- x[itr, ] # train
ytr <- y[itr]
xte <- x[-itr, ] # test
yte <- y[-itr]

# take prototype candidates identical to training points:
z <- xtr
dxz <- dist2(xtr, z)
# run protoclass:
prot <- protoclass(dxz=dxz, y=ytr, eps=2, lambda=1/n)
\dontrun{
plot(prot,xtr,y=1+(ytr=="A"))
}
# get predictions on test data:
pred1 <- predict(prot, xte, z=xtr)
# get predictions on test data using pairwise distances:
pred2 <- predictwithd.protoclass(prot, dist2(xte, z))
}
\seealso{
  \code{\link{predict.protoclass}}
}

