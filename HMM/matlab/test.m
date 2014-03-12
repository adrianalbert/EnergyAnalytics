% TEST test script for HHM package

% References: Hidden Markov Models by Cappe, Moulines, Rydden
% Springer series in statistics

% by Aurelien Garivier, CNRS & Telecom ParisTech
% last revision January 30, 2012


n = 1000;
nu = [0, 1];
Q = [0.8, 0.2; 0.1, 0.9];
g = [0.25 0.25 0.25 0.25; 0.05 0.05 0.45 0.45];

tic; [x,y] = HMMsample(nu, Q, g, n); toc
tic; [phi, c] = HMMfilter(y, nu, Q, g); toc
[~, xh] = max(phi, [], 1);
mean(x==xh) % estimation performance of the filter

tic; beta = HMMsmoother(y, Q, g, c); toc
psi = phi .* beta; % posterior probabilities
[~, xh] = max(psi, [], 1); % MAP
mean(x==xh) % estimation performance of the marginal MAP

tic; [Qh, gh] = HMMbaumwelch(y, nu); toc
Q, Qh
g, gh

%compare real and maximum log-likelihood
[~, ch] = HMMfilter(y, nu, Qh, gh);
sum(log(ch))-sum(log(c))