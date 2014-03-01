function beta = HMMsmoother(y, Q, g, c)
% HMMsmoother compute smoothing coefficients
%
%  in : y = vector of observations of size n, with values between 1 and r
%       Q = transition matrix of size k
%       g = emission matrix of size k x r
%       c = vector computed by function filter: c(t) = P(Y(t) = y(t)| Y(1:t-1)=y(1:t-1))
% out : beta = smoothing factors: 
%         beta(x,t) = P(Y(t+1:n)=y(t+1:n) | X(t)=x) / P(Y(t+1:n)=y(t+1:n) | Y(1:t)=y(1:t)) for 1<=x<=k and 1<=t<=n
%       permits to compute the posterior distribution of the hidden states as post = phi .* beta:
%         post(x,t) = P(X(t)=x | Y(1:n)=y(1:n))
%
% Example : 
%       n = 100;
%       nu = [0, 1];
%       Q = [0.8, 0.2; 0.1, 0.9];
%       g = [0.25 0.25 0.25 0.25; 0.05 0.05 0.45 0.45];
%       [x,y] = HMMsample(nu, Q, g, n);
%       [phi, c] = HMMfilter(y, nu, Q, g);
%       beta = HMMsmoother(y, Q, g, c);
%       psi = phi .* beta; % posterior probabilities
%       [~, xh] = max(psi, [], 1); % MAP
%       mean(x==xh) % performance of the filter

% References: Hidden Markov Models by Cappe, Moulines, Rydden
% Springer series in statistics

% by Aurelien Garivier, CNRS & Telecom ParisTech
% last revision January 30, 2012


n = length(y);
beta = ones(size(Q, 1), n);
for t=(n-1):-1:1
  beta(:, t) = Q * (g(:, y(t+1)) .* beta(:, t+1)) / c(t+1);
end
