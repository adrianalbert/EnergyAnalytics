// written by Aurelien Garivier, CNRS & Telecom Paristech
// January 2012
//
// Baum-Welch algorithm for discrete Hidden Markov models
// see http://www.telecom-paristech.fr/~garivier/code/index.html

function [x,y] = HMMsample(nu, Q, g, n)
// HMMsample sample a trajectory from a hidden markov chain
//
//  in : nu = initial distribution as vector of size k
//       Q = transition matrix of size k
//       n = positive integer
// out :  [x,y] = sample trajectory of size n of a HMM defined by (nu, Q, g):
//       x = sample trajectory of size n of a Markov Chain with initial distribution nu and transition matrix Q
//       y = observations such that the conditionnal distribution of y(k)
//       given x(k) is g(x(k), :)
//
// Example : 
//       n = 100;
//       nu = [0, 1];
//       Q = [0.8, 0.2; 0.1, 0.9];
//       g = [0.25 0.25 0.25 0.25; 0.05 0.05 0.45 0.45];
//       [x,y] = HMMsample(nu, Q, g, n);

  cQ = cumsum(Q, 2); cg = cumsum(g, 2);
  x = zeros(1, n); y = zeros(1, n);
  x(1) = 1+sum(rand()>cumsum(nu));
  y(1) = 1+sum(rand()>cg(x(1), :));

  for j=2:n
      x(j) = 1+sum(rand()>cQ(x(j-1), :));
      y(j) = 1+sum(rand()>cg(x(j), :));
  end
endfunction


function [phi,c] = HMMfilter(y,nu,Q,g)
// HMMfilter filter of the observation sequence given HMM parameters
//
//  in : y = vector of observations of size n, with values between 1 and r
//       nu = initial distribution as vector of size k
//       Q = transition matrix of size k
//       g = emission matrix of size k x r
// out : phi = filter(x,t) = P(X(t)=x | Y(1:t)=y(1:t)) for 1<=x<=k and 1<=t<=n
//       c(t) = P(Y(t) = y(t)| Y(1:t-1)=y(1:t-1))
//
// Example : 
//       n = 100;
//       nu = [0, 1];
//       Q = [0.8, 0.2; 0.1, 0.9];
//       g = [0.25 0.25 0.25 0.25; 0.05 0.05 0.45 0.45];
//       [x,y] = HMMsample(nu, Q, g, 20);
//       [phi, c] = HMMfilter(y, nu, Q, g);
//       [~, xh] = max(phi, [], 1); // pseudo-MAP
//       mean(x==xh) // performance of the filter

  n = max(size(y));
  phi = zeros(size(Q,1),n);
  c = zeros(1,n);

  Z = nu' .* g(:,y(1));
  c(1) = sum(Z);
  phi(:,1) = Z/c(1);

  for t = 2:n
    Z = (phi(:,t-1)' * Q)' .* g(:,y(t));
    c(t) = sum(Z);
    phi(:,t) = Z/c(t);
  end;
endfunction


function [betaa] = HMMsmoother(y,Q,g,c)
// HMMsmoother compute smoothing coefficients
//
//  in : y = vector of observations of size n, with values between 1 and r
//       Q = transition matrix of size k
//       g = emission matrix of size k x r
//       c = vector computed by function filter: c(t) = P(Y(t) = y(t)| Y(1:t-1)=y(1:t-1))
// out : beta = smoothing factors: 
//         betaa(x,t) = P(Y(t+1:n)=y(t+1:n) | X(t)=x) / P(Y(t+1:n)=y(t+1:n) | Y(1:t)=y(1:t)) for 1<=x<=k and 1<=t<=n
//       permits to compute the posterior distribution of the hidden states as post = phi .* beta:
//         post(x,t) = P(X(t)=x | Y(1:n)=y(1:n))
//
// Example : 
//       n = 100;
//       nu = [0, 1];
//       Q = [0.8, 0.2; 0.1, 0.9];
//       g = [0.25 0.25 0.25 0.25; 0.05 0.05 0.45 0.45];
//       [x,y] = HMMsample(nu, Q, g, n);
//       [phi, c] = HMMfilter(y, nu, Q, g);
//       betaa = HMMsmoother(y, Q, g, c);
//       psi = phi .* beta; // posterior probabilities
//       [~, xh] = max(psi, [], 1); // MAP
//       mean(x==xh) // performance of the filter

  n = max(size(y));
  betaa = ones(size(Q,1),n);
  for t = n-1:-1:1
    betaa(:,t) = Q * (g(:, y(t+1)) .* betaa(:, t+1)) / c(t+1);
  end;
endfunction


function [Q,g,l] = HMMbaumwelch(y,nu,tol,maxIt)
// HMMbaumwelch compute maximum likehood estimate using Expectation-Maximization
// iterations
//
//  in : y = vector of observations 
//       nu = initial distribution of the hidden chain
//       tol = tolerance for the stopping criterion
//       maxIt = maximal number of iterations
// out : Q = estimate of the transition matrix of the hidden markov process
//       g = estimated probabilities of transition: gh(x,y) = estimate of P(Y=y | X=x) for 1<=x<=k
//       l = log-likelihood of y for parameters Q and g
//
// Example :
//       n = 10000;
//       nu = [0, 1];
//       Q = [0.8, 0.2; 0.1, 0.9];
//       g = [0.25 0.25 0.25 0.25; 0.05 0.05 0.45 0.45];
//       [x,y] = HMMsample(nu, Q, g, n);
//       [Qh, gh] = HMMem(y, nu);
//       // compare estimates with truth: note that the order of the hidden
//       // states may not be preserved
//       Q, Qh
//       g, gh

  [lhs,rhs]=argn(0)
  if rhs<4 then maxIt = 100; end;
  if rhs<3 then tol = 0.0001; end;

  k = max(size(nu)); r = max(y); n = max(size(y));
  Y = zeros(n, r); Y(sub2ind([n, r], 1:n, y))=1;
  Q = rand(k,k); Q = Q ./(sum(Q,2)*ones(1,k)); // random initial transition matrix
  g = rand(k, r); g = g ./ (sum(g, 2)*ones(1, r));
  it = 0; oldQ = Q; oldg = g+tol+1;
  while ((norm(oldQ(:)-Q(:), 1) + norm(oldg-g, 1) > tol) & (it<maxIt))
   it = it+1;
   // compute the posterior distribution for the current parameters
   [phi,c] = HMMfilter(y,nu,Q,g);
   betaa = HMMsmoother(y,Q,g,c);
   post = phi .*betaa;

   // expectation of the number of transitions under the current parameters
   N = Q .*(phi(:, 1:$-1    )*(betaa(:,2:$   ) .*g(:, y(2:$)  )./ (ones(k,1)*c(2:$   )))'); 
   // expectation of the number of emissions
   M = post*Y;

   // re-estimation
   oldQ = Q; oldg = g;
   Q = N ./(sum(N,2)*ones(1,k));
   g = M ./(sum(M,2)*ones(1,r));
  end;
  l = sum(log(c));
endfunction

function [] = benchmark(nbReps, scenario, finish)
// benchmark: BENCHMARK for the comparison with other programming languages
// see http://www.telecom-paristech.fr/~garivier/code/index.html

  [lhs,rhs]=argn(0)
  if rhs<3 then finish = 0; end
  if rhs<2 then scenario = 1; end
  if rhs<1 then nbReps = 10000; end
  
  if (scenario == 1)
    nu = [0,1];
    Q = [0.8,0.2;0.1,0.9];
    g = [0.25,0.25,0.25,0.25;0.1,0.1,0.4,0.4];
  else
    k = 50; e = 0.1;
    nu = zeros(1, k); nu(1)=1;
    Q = (1-2*e)*eye(k,k) + e*diag([2, ones(1, k-2)], 1) + e*diag([ones(1, k-2), 2], -1);
    g = Q;
  end
  disp(["running " , string(nbReps), " repetitions of scenario ", string(scenario)])

  n = 10000;
  tps = zeros(nbReps,3);
  for rep = 1:nbReps, disp(rep)
     tic;
     [x,y] = HMMsample(nu,Q,g,n);
     tps(rep,1) = toc();

     // compute filter and smoother
     tic;
     [phi,c] = HMMfilter(y,nu,Q,g);
     betaa = HMMsmoother(y,Q,g,c);
     post = phi .* betaa;
     tps(rep,2) = toc();
     l0 = sum(log(c)); // log-likelihood of the observations

     // estimate parameters from observations only
     nbTrials = 10; // number of searchs for a local maximum
     bestl = -1e300;
     tic;
     for j = 1:nbTrials
        [Qn,gn,l] = HMMbaumwelch(y,nu);
        if l>bestl then
          Qh = Qn;  gh = gn;  bestl = l;
        end;
    end;
    tps(rep,3) = toc();
    bestl, Qh, gh

  // compare likelihoods between estimates and real parameters
  //[bestl, l0]
  
  end;
  disp(tps)
  disp(mean(tps,1))
  if finish then quit; end
endfunction

