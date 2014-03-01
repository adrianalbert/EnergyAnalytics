function [phi, c] = HMMfilter(y, nu, Q, g)
% HMMfilter filter of the observation sequence given HMM parameters
%
%  in : y = vector of observations of size n, with values between 1 and r
%       nu = initial distribution as vector of size k
%       Q = transition matrix of size k
%       g = emission matrix of size k x r
% out : phi = filter(x,t) = P(X(t)=x | Y(1:t)=y(1:t)) for 1<=x<=k and 1<=t<=n
%       c(t) = P(Y(t) = y(t)| Y(1:t-1)=y(1:t-1))
%
% Example : 
%       n = 100;
%       nu = [0, 1];
%       Q = [0.8, 0.2; 0.1, 0.9];
%       g = [0.25 0.25 0.25 0.25; 0.05 0.05 0.45 0.45];
%       [x,y] = HMMsample(nu, Q, g, 20);
%       [phi, c] = HMMfilter(y, nu, Q, g);
%       [~, xh] = max(phi, [], 1); % pseudo-MAP
%       mean(x==xh) % performance of the filter

% References: Hidden Markov Models by Cappe, Moulines, Rydden
% Springer series in statistics

% by Aurelien Garivier, CNRS & Telecom ParisTech
% last revision January 30, 2012


n = length(y);
phi = zeros(size(Q, 1), n);
c = zeros(1, n);

Z = nu'.*g(:, y(1));
c(1)= sum(Z);
phi(:, 1) = Z/c(1);

for t=2:n
  Z = (phi(:, t-1)' * Q)' .* g(:, y(t));
  c(t) = sum(Z);
  phi(:, t) = Z / c(t);
end