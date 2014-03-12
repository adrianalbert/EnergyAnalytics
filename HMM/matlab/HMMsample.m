function [x,y] = HMMsample(nu, Q, g, n)
% HMMsample sample a trajectory from a hidden markov chain
%
%  in : nu = initial distribution as vector of size k
%       Q = transition matrix of size k
%       n = positive integer
% out :  [x,y] = sample trajectory of size n of a HMM defined by (nu, Q, g):
%       x = sample trajectory of size n of a Markov Chain with initial distribution nu and transition matrix Q
%       y = observations such that the conditionnal distribution of y(k)
%       given x(k) is g(x(k), :)
%
% Example : 
%       n = 100;
%       nu = [0, 1];
%       Q = [0.8, 0.2; 0.1, 0.9];
%       g = [0.25 0.25 0.25 0.25; 0.05 0.05 0.45 0.45];
%       [x,y] = HMMsample(nu, Q, g, n);

% References: Hidden Markov Models by Cappe, Moulines, Rydden
% Springer series in statistics

% by Aurelien Garivier, CNRS & Telecom ParisTech
% last revision January 30, 2012

cQ = cumsum(Q, 2); cg = cumsum(g, 2);
x = zeros(1, n); y = zeros(1, n);
x(1) = 1+sum(rand>cumsum(nu));
y(1) = 1+sum(rand>cg(x(1), :));

for j=2:n
    x(j) = 1+sum(rand>cQ(x(j-1), :));
    y(j) = 1+sum(rand>cg(x(j), :));
end