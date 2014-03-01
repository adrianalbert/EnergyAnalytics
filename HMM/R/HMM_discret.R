## Created on Jan 26, 2012
##
## @author: Aurelien Garivier, CNRS and Telecom ParisTech
## 
## Baum-Welch algorithm for discrete Hidden Markov models
## 
## @reference: Hidden Markov Models by Cappe, Moulines, Rydden
## Springer series in statistics
## see also http://www.telecom-paristech.fr/~garivier/code/index.html


## sample a trajectory from a hidden markov chain
## in:   nu = initial distribution as vector of size k
##       Q = transition matrix of size k
##       n = positive integer
## out:  (x,y) = sample trajectory of size n of a HMM defined by (nu, Q, g):
##       x = sample trajectory of size n of a Markov Chain with initial distribution nu and transition matrix Q
##       y = observations such that the conditionnal distribution of y[k]
##       given x[k] is g(x[k], :)
HMMsample <- function(nu, Q, g, n){
  cQ <- t(apply(Q, 1, cumsum));
  cg <- t(apply(g, 1, cumsum));
  x <- array(0, n);
  y <- array(0, n);
  x[1] <- 1+sum(as.numeric(runif(1)>cumsum(nu)));
  y[1] <- 1+sum(as.numeric(runif(1)>cg[x[1],]));
  for (t in 2:n){
    x[t] <- 1+sum(as.numeric(runif(1)>cQ[x[t-1],])); 
    y[t] <- 1+sum(as.numeric(runif(1)>cg[x[t],]));
    #remark: it is slightly simpler, but also slightly slower to use:
    # x[t] <- sample(1:k, 1, prob = Q[x[t-1],])
    # y[t] <- sample(1:r, 1, prob = g[x[t],]);
    
  }
  list(x=x,y=y);
}

## HMM filtering of an observation sequence, given hmm parameters
## in:   y = vector of observations, assumed to be in range(g.shape[1])
##       nu = initial distribution as vector of size k
##       Q = transition matrix of size k x k
##       g = emission matrix with k rows
## out:  phi = filter: P(x[t]=x | y[1:t]=y[1:t]) for 1<=x<k and 0<=t<n
##       c(t) = conditional likelihood: P(Y[t] = y[t]| Y[1:t-1]=y[1:t-1])
HMMfilter <- function(y, nu, Q, g){
  n <- length(y);
  dims <- dim(Q);
  res <- list();
  res$phi <- matrix(nrow = dims[1], ncol = n);
  res$c <- array(0, n);

  Z <- nu*g[,y[1]];
  res$c[1] <- sum(Z);
  res$phi[,1] <- Z/res$c[1];

  for (t in 2:n){
    Z <- (res$phi[,t-1] %*% Q) * g[, y[t]]
    res$c[t] <- sum(Z);
    res$phi[,t] <- Z / res$c[t];
  }
  res;
}

## HMM filtering of an observation sequence, given hmm parameters
## in:   y = vector of observations, assumed to be in range(Q.shape[0])
##       Q = transition matrix of size k x k
##       g = emission matrix with k rows
##       c = conditional likelihoods, computed by HMMfilter
## out:  beta = smoothing factors: P(y[t+1:n]=y[t+1:n] | X[t]=x) / P(Y[t+1:n]=y[t+1:n] | Y[0:t]=y[1:t]) for 0<=x<k and 1<=t<n
##       permits to compute the posterior distribution of the hidden states 
##       P(X[t]=x | Y[0:n]=y[0:n])  as post = phi * beta
HMMsmoother <- function(y, Q, g, c){
  n<- length(y);
  dims <- dim(Q);
  beta <- matrix(1, nrow=dims[1], ncol=n);
  for (t in seq(n-1, 1, -1)){
    beta[,t] = Q %*% (g[,y[t+1]] * beta[, t+1]) / c[t+1];
  }
  beta;
}


# to load the C filter/smoother, you must first compile them using "make"
dyn.load("HMM.so")

# R wrapper for the C filter calls
HMMfilter_C = function(y, nu, Q, g){
  n <- length(y);
  N <- dim(Q)[1];
  M <- dim(g)[2];

  #res <- .C("HMMfilter", as.integer(N), as.integer(M), as.integer(n), as.double(y), as.double(nu), as.double(Q), as.double(g), as.double(matrix(0, nrow = N, ncol = n)), as.double(array(0, n)));
  res <- .C("HMMfilter", as.integer(N), as.integer(M), as.integer(n), as.double(y), as.double(nu), as.double(Q), as.double(g), as.double(array(0, N*n)), as.double(array(0, n)));
  list(phi = matrix(res[[8]], nrow = N, ncol = n), c = res[[9]]);
}

# R wrapper for the C smoother calls
HMMsmoother_C = function(y, Q, g, c){
  n <- length(y);
  N <- dim(Q)[1];
  M <- dim(g)[2];

  res <- .C("HMMsmoother", as.integer(N), as.integer(M), as.integer(n), as.double(y), as.double(Q), as.double(g), as.double(c), as.double(matrix(0, nrow = N, ncol = n)))[[8]];
  matrix(res, nrow = N, ncol = n);
}

## compute maximum likehood estimate using Expectation-Maximization iterations
## in:   y = vector of observations 
##       nu = initial distribution of the hidden chain
##       tol = tolerance for the stopping criterion
##       maxIt = maximal number of iterations
## out:  Q = estimate of the transition matrix of the hidden markov process
##       g = estimated probabilities of transition: g(x,y) = estimate of P(Y=y | X=x) for 0<=x<k
##       l = log-likelihood of y for parameters Q and g  
HMMbaumwelch <- function(y, nu, tol = 1e-4, maxIt = 100){
# requires global variable myfilter to be set, either to HMMfilter or to HMMfilter_C
# idem for mysmoother
  k <- length(nu);
  r <- max(y);
  n <- length(y);
  Y <- matrix(0, nrow=n, ncol=r);   for (t in 1:n){ Y[t, y[t]]=1; }
  Q <- matrix(runif(k*k), nrow=k);  Q <- Q / apply(Q, 1, sum);
  g <- matrix(runif(k*r), nrow=k);  g <- g / apply(g, 1, sum);
  it <- 0; oldQ <- Q-tol; oldg <- g-tol;
  while ((sum(abs((oldQ-Q))) + sum(abs(oldg-g)) > tol) & (it<maxIt)){
    it <- it+1;
    f <- myfilter(y, nu, Q, g);
    beta <- mysmoother(y, Q, g, f$c);
    post <- f$phi * beta;
  
    N <- Q * (f$phi[,1:(n-1)] %*% t(beta[, 2:n] * g[, y[2:n]]  / (matrix(1, nrow=k, ncol=1)%*%f$c[2:n])));
    M <- post %*% Y;
    
    oldQ <- Q; oldg <- g;
    Q <- N / apply(N, 1, sum);
    g <- M / apply(M, 1, sum);
  }
  
  res <- list();
  res$Q <- Q;
  res$g <- g;
  res$l <- sum(log(f$c));
  res$it <- it;
  res;
}

## BENCHMARK for the comparison with other programming languages
## see http://www.telecom-paristech.fr/~garivier/code/index.html  
benchmark <- function(nbReps = 10, scenario = 1, withC = 0){
	if (scenario==1){
		nu <- c(0, 1);
		Q <- matrix(c(0.8, 0.2, 0.1, 0.9), byrow=TRUE, nrow=2);
		g <- matrix(c(0.25, 0.25, 0.25, 0.25, 0.1, 0.1, 0.4, 0.4), byrow=TRUE, nrow=2);
	}
	else{
		k <- 50; e <- 0.1;
		nu <- c(1, rep(0, k-1));
		Q <- matrix(0, nrow=k, ncol=k);
		for(i in 1:k) Q[i,i] = 1-2*e;
		for(i in 2:(k-1)){ Q[i, i-1] = e; Q[i, i+1] = e;}
		Q[1, 2] = 2*e; Q[k, k-1] = 2*e;
		g <- Q;	
	}
	if (withC){
		myfilter <<- HMMfilter_C;
		mysmoother <<- HMMsmoother_C;
	}
	else{
		myfilter <<- HMMfilter;
		mysmoother <<- HMMsmoother;
	}
	wit = c("with pure R filtering/smoothing", "with C filtering/smoothing");
	print(paste("running",nbReps,"repetitions of scenario",scenario,wit[1+(withC>0)], sep=" "));
	
	n <- 10000;	
	tps <- matrix(0, nrow = nbReps, ncol=3);
	for(rep in 1:nbReps){ print(rep);
		pre.time <- Sys.time()
		traj <- HMMsample(nu, Q, g, n)
		tps[rep, 1] <- Sys.time()-pre.time

		pre.time <- Sys.time()
		f <- myfilter(traj$y, nu, Q, g)
		beta <- mysmoother(traj$y, Q, g, f$c)
		psi <- f$phi * beta;
		tps[rep, 2] <- Sys.time()-pre.time

		nbTrials <- 10;
		pre.time <- Sys.time()
		bestl = -Inf;
		for(j in 1:nbTrials){
			estim <- HMMbaumwelch(traj$y, nu)
			if (estim$l > bestl){
				Qh <- estim$Q
				gh <- estim$g
				bestl <- estim$l
			}
		}
		tps[rep, 3] <-Sys.time()-pre.time
	}
	print(tps)
	print(colMeans(tps))
}

benchmarkPackageHMM <- function(nbReps = 10, scenario = 1){
	require(HMM)
	if (scenario==1){
		nu <- c(0, 1);
		Q <- matrix(c(0.8, 0.2, 0.1, 0.9), byrow=TRUE, nrow=2);
		g <- matrix(c(0.25, 0.25, 0.25, 0.25, 0.1, 0.1, 0.4, 0.4), byrow=TRUE, nrow=2);
	}
	else{
		k <- 50; e <- 0.1;
		nu <- c(1, rep(0, k-1));
		Q <- matrix(0, nrow=k, ncol=k);
		for(i in 1:k) Q[i,i] = 1-2*e;
		for(i in 2:(k-1)){ Q[i, i-1] = e; Q[i, i+1] = e;}
		Q[1, 2] = 2*e; Q[k, k-1] = 2*e;
		g <- Q;	
	}
	hmm = initHMM(1:dim(Q)[1], 1:dim(g)[2], transProbs=Q, emissionProbs=g)
	print(paste("running",nbReps,"repetitions of scenario",scenario,"using package HMM", sep=" "));
	
	n <- 10000;	
	tps <- matrix(0, nrow = nbReps, ncol=3);
	for(rep in 1:nbReps){ print(rep);
		pre.time <- Sys.time()
		traj <- simHMM(hmm, n)
		tps[rep, 1] <- Sys.time()-pre.time

		pre.time <- Sys.time()
		psi <- posterior(hmm, traj$observation)
		tps[rep, 2] <- Sys.time()-pre.time

		nbTrials <- 10;
		pre.time <- Sys.time()
		bestl = -Inf;
		for(j in 1:nbTrials){
			estim <- baumWelch(hmm, traj$observation, maxIterations=100, delta=1E-9, pseudoCount=0)
		}
		tps[rep, 3] <-Sys.time()-pre.time
	}
	print(tps)
	print(colMeans(tps))
}

