'''
Created on Jan 26, 2012

@author: Aurelien Garivier, CNRS and Telecom ParisTech

Baum-Welch algorithm for discrete Hidden Markov models

@reference: Hidden Markov Models by Cappe, Moulines, Rydden
Springer series in statistics
see also http://www.telecom-paristech.fr/~garivier/code/index.html
'''

from random import random, choice
import numpy as np
import sys
import time
import HMM_C # needs the C functions to be compiled: use "make"

def randm(p):
  '''
sample from discrete distribution
in:   p = vector of probabilities, assumed to sum to 1
  '''
  res = 0
  q = p[0]
  u = random()
  while(u>q):
    res += 1
    q += p[res]
  return res

def HMMsample(nu, Q, g, n):
  '''
sample a trajectory from a hidden markov chain
in:   nu = initial distribution as vector of size k
      Q = transition matrix of size k
      n = positive integer
out:  (x,y) = sample trajectory of size n of a HMM defined by (nu, Q, g):
      x = sample trajectory of size n of a Markov Chain with initial distribution nu and transition matrix Q
      y = observations such that the conditionnal distribution of y[k]
      given x[k] is g(x[k], :)
  '''
  x = np.zeros(n)
  y = np.zeros(n)
  x[0] = randm(nu)
  y[0] = randm(g[x[0], :])
  for j in range(1,n):
    x[j] = randm(Q[x[j-1], :])
    y[j] = randm(g[x[j], :])
  return (x,y)

def HMMfilter(y, nu, Q, g):
  '''
HMM filtering of an observation sequence, given hmm parameters
in:   y = vector of observations, assumed to be in range(g.shape[1])
      nu = initial distribution as vector of size k
      Q = transition matrix of size k x k
      g = emission matrix with k rows
out:  phi = filter: P(x[t]=x | y[0:t]=y[0:t]) for 0<=x<k and 0<=t<n
      c(t) = conditional likelihood: P(Y[t] = y[t]| Y[0:t-1]=y[0:t-1])
  '''  
  n = y.size
  k = Q.shape[0]
  phi = np.zeros([k, n])
  c = np.zeros(n)
  
  Z = nu * g[:, y[0]]
  c[0] = sum(Z)
  phi[:,0] = Z/c[0]
  
  for t in range(1,n):
    Z = np.dot(phi[:, t-1], Q) * g[:, y[t]];
    c[t] = sum(Z)
    phi[:, t] = Z / c[t]
  
  return (phi,c)


def HMMsmoother(y, Q, g, c):
  '''
HMM filtering of an observation sequence, given hmm parameters
in:   y = vector of observations, assumed to be in range(Q.shape[0])
      Q = transition matrix of size k x k
      g = emission matrix with k rows
      c = conditional likelihoods, computed by HMMfilter
out:  beta = smoothing factors: P(y[t+1:n]=y[t+1:n] | X[t]=x) / P(Y[t+1:n]=y[t+1:n] | Y[0:t]=y[1:t]) for 0<=x<k and 1<=t<n
      permits to compute the posterior distribution of the hidden states 
      P(X[t]=x | Y[0:n]=y[0:n])  as post = phi * beta
  '''
  n = y.size
  beta = np.ones([Q.shape[0], n])
  for t in range(n-2, -1, -1):
    beta[:,t] = np.dot(Q, g[:, y[t+1]] * beta[:, t+1]) / c[t+1]
  return beta

def HMMbaumwelch(y, nu, tol=1e-4, maxIt = 100):
  '''
compute maximum likehood estimate using Expectation-Maximization
iterations
in:   y = vector of observations 
      nu = initial distribution of the hidden chain
      tol = tolerance for the stopping criterion
      maxIt = maximal number of iterations
out:  Q = estimate of the transition matrix of the hidden markov process
      g = estimated probabilities of transition: g(x,y) = estimate of P(Y=y | X=x) for 0<=x<k
      l = log-likelihood of y for parameters Q and g  
  '''  
  global myfilter, mysmoother # should be either HMMfilter/HMMsmoother, or HMM_C.HMMfilter/HMM_C.HMMsmoother
  k = nu.size; r = 1+max(y); n = y.size
  Y = np.zeros([n, r]); 
  Y[range(n), np.int_(y)] = 1;
  Q = np.random.rand(k, k)
  g = np.random.rand(k, r)
  for j in range(k):
    Q[j,:] = Q[j,:] / sum(Q[j,:])
    g[j,:] = g[j,:] / sum(g[j,:])
  it = 0; oldQ = Q; oldg = g + tol +1
  
  while (sum(sum(abs(oldQ[:]-Q[:]))) + sum(sum(abs(oldg-g))) > tol) & (it<maxIt):
    it+=1
    (phi, c) = myfilter(y, nu, Q, g)
    beta = mysmoother(y, Q, g, c)
    post = phi * beta
    N = Q * (np.dot(phi[:, 0:-1], np.transpose(beta[:, 1:]*g[:, np.int_(y[1:])]/np.tile(c[1:], [k, 1]))))
    M = np.dot(post, Y)

    oldQ = Q.copy(); oldg = g.copy()
    for j in range(k):
      Q[j,:] = N[j,:] / sum(N[j,:])
      g[j,:] = M[j,:] / sum(M[j,:])
  l = sum(np.log(c))
  return (Q, g, l)

def benchmark(nbReps=10, scenario=1, withC=0):
  '''
BENCHMARK for the comparison with other programming languages
see http://www.telecom-paristech.fr/~garivier/code/index.html  
  '''  
  global myfilter, mysmoother # should be either HMMfilter/HMMsmoother, or HMM_C.HMMfilter/HMM_C.HMMsmoother
  
  if scenario == 1:
    nu = np.array([0., 1.])
    Q = np.array([0.8, 0.2, 0.1, 0.9])
    Q.shape = [2,2]
    g = np.array([0.25, 0.25, 0.25, 0.25, 0.1, 0.1, 0.4, 0.4])
    g.shape = [2, 4]
  else:
    k = 50
    e = 0.1;
    nu = np.zeros(k); nu[0] = 1.;
    Q = np.zeros([k, k]);
    for i in range(k):
      Q[i, i] = 1 - 2*e
      if i>0:
        if i<k-1:
          Q[i, i-1] = e
          Q[i, i+1] = e
        else:
          Q[i,i-1] = 2*e
      else:
        Q[i, i+1] = 2*e
    g = Q;
 
  # choose between pure python and C filtering/smoothing
  if withC:
    myfilter = HMM_C.HMMfilter
    mysmoother = HMM_C.HMMsmoother
  else:
    myfilter = HMMfilter
    mysmoother = HMMsmoother
  
  n = 100
  tps = np.zeros([nbReps, 3])
      
  for rep in range(nbReps):
    print rep;
    # create the sample
    t0 = time.time()
    (x, y) = HMMsample(nu, Q, g, n)
    tps[rep, 0] = time.time() - t0;
    
    # filtering / smoothing
    t0 = time.time()
    (phi, c) = myfilter(y, nu, Q, g)
    beta = mysmoother(y, Q, g, c)
    post = phi * beta;
    tps[rep,1] = time.time()-t0;
      
    # estimate parameters from observations only 
    nbTrials = 10; # number of searchs for a local maximum
    bestl=-1e300;
    t0=time.time() 
    for j in range(nbTrials):
      (Qn, gn, l) = HMMbaumwelch(y, nu)
      if l>bestl:
        Qh = Qn; gh = gn; bestl = l;
    tps[rep, 2] = time.time() - t0;
  print(tps)
  print(np.mean(tps, 0))  

if __name__ == "__main__":
  if len(sys.argv)<4:
    withC = 0
  else:
    withC = int(sys.argv[3]);
  if len(sys.argv)<3:
    scenario = 1
  else:
    scenario = int(sys.argv[2]);    
  if len(sys.argv)<2:
    nbReps = 10
  else:
    nbReps = int(sys.argv[1]); 
  wit = [" in pure python", " with C filtering/smoothing"]
  print "running "+str(nbReps)+ " repetitions of scenario "+ str(scenario) + wit[withC>0] + "..."
  benchmark(nbReps, scenario, withC)

  
