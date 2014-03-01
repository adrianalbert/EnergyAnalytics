% BENCHMARK for the comparison with other programming languages
% see http://www.telecom-paristech.fr/~garivier/code/index.html
%
% by Aurelien Garivier, CNRS & Telecom ParisTech
% last revision January 30, 2012

nbReps = 2;
tps = zeros(nbReps, 3);

nu = [1, 0];
Q = [0.8, 0.2; 0.1, 0.9];
g = [0.25 0.25 0.25 0.25; 0.1 0.1 0.4 0.4];
n = 10000;

for rep = 1:nbReps, rep
% create the sample
    tic; 
    [x,y] = hmmgenerate(n, Q, g); % assumes x(0) == 1
    tps(rep, 1) = toc;

% filtering / smoothing    
    tic; 
    post = hmmdecode(y, Q, g); % assumes x(0) == 1
    tps(rep, 2) = toc;

% estimate parameters from observations only
	k = size(Q, 1); r = size(g, 2);
    nbTrials = 10; 
    bestl=-inf;
    tic; 
    for j=1:nbTrials
      Q0 = rand(k); Q = Q ./ (sum(Q, 2)*ones(1, k));
      g0 = rand(k, r); g = g ./ (sum(g, 2)*ones(1, r));
      [Qn, gn] = hmmtrain(y, Q0, g0, 'Maxiterations', 100); 
      % should check if estimate is better, using likelihood...
    end
    tps(rep, 3) = toc;
end
disp(tps)
disp(mean(tps, 1))
exit
