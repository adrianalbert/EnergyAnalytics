function [Q, g, l] = HMMbaumwelch(y, nu, tol, maxIt, Q, g)
% HMMbaumwelch compute maximum likehood estimate using Expectation-Maximization
% iterations
%
%  in : y = vector of observations 
%       nu = initial distribution of the hidden chain
%       tol = tolerance for the stopping criterion
%       maxIt = maximal number of iterations
% out : Q = estimate of the transition matrix of the hidden markov process
%       g = estimated probabilities of transition: gh(x,y) = estimate of P(Y=y | X=x) for 1<=x<=k
%       l = log-likelihood of y for parameters Q and g
%
% Example :
%       n = 10000;
%       nu = [0, 1];
%       Q = [0.8, 0.2; 0.1, 0.9];
%       g = [0.25 0.25 0.25 0.25; 0.05 0.05 0.45 0.45];
%       [x,y] = HMMsample(nu, Q, g, n);
%       [Qh, gh] = HMMem(y, nu);
%       % compare estimates with truth: note that the order of the hidden
%       % states may not be preserved
%       Q, Qh
%       g, gh

% References: Hidden Markov Models by Cappe, Moulines, Rydden
% Springer series in statistics

% by Aurelien Garivier, CNRS & Telecom ParisTech
% last revision February 7, 2012

global myfilter mysmoother % should be either HMMfilter/HMMsmoother, or HMMfilter_C/HMMsmoother_C

if nargin<4, maxIt = 100; end
if nargin<3, tol = 1e-4; end

k = length(nu); r = max(y); n = length(y);
Y = zeros(n, r); Y(sub2ind([n, r], 1:n, y))=1;

% if they are not provided, sample random initial transition and emission matrices
if nargin<5, Q = rand(k); Q = Q ./ (sum(Q, 2)*ones(1, k)); end 
if nargin<6, g = rand(k, r); g = g ./ (sum(g, 2)*ones(1, r)); end 

it = 0; oldQ = Q; oldg = g+tol+1;
while ((norm(oldQ(:)-Q(:), 1) + norm(oldg-g, 1) > tol) && (it<maxIt))
  it = it + 1;
  % compute the posterior distribution for the current parameters
  [phi, c] = myfilter(y, nu, Q, g);
  beta = mysmoother(y, Q, g, c);
  post = phi .* beta;
  
  % expectation of the number of transitions under the current parameters
  N =Q.*(phi(:, 1:(end-1))*(beta(:, 2:end).*g(:, y(2:end))./(ones(k, 1)*c(2:end)))');   
  % expectation of the number of emissions
  M = post * Y;

  % re-estimation
  oldQ = Q; oldg = g;
  Q = N ./ (sum(N, 2) * ones(1, k));
  g = M ./ (sum(M, 2) * ones(1, r));
end
l = sum(log(c));
