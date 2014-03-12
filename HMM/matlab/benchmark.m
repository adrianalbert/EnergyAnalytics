% BENCHMARK for the comparison with other programming languages
% see http://www.telecom-paristech.fr/~garivier/code/index.html
%
% by Aurelien Garivier, CNRS & Telecom ParisTech
% last revision January 30, 2012


function benchmark(nbReps, scenario, withC, finish)
	global myfilter mysmoother

	if nargin<1, nbReps = 10; end
	if nargin<2, scenario = 1; end
	if nargin<3, withC = 0; end
	if nargin<4, finish = 0; end
	wit = {' with pure matlab filtering/smoothing', ' with C filtering/smoothing'};
	disp(['running ', num2str(nbReps), ' repetitions of scenario ', num2str(scenario),  wit{withC+1}]);

	if scenario==1
		nu = [0, 1];
		Q = [0.8, 0.2; 0.1, 0.9];
		g = [0.25 0.25 0.25 0.25; 0.1 0.1 0.4 0.4];
	else
		k = 50; e = 0.1;
		nu = zeros(1, k); nu(1)=1;
		Q = (1-2*e)*eye(k) + e*diag([2, ones(1, k-2)], 1) + e*diag([ones(1, k-2), 2], -1);
		g = Q;
	end

	if withC
		myfilter = @HMMfilter_C;
		mysmoother = @HMMsmoother_C;
	else
		myfilter = @HMMfilter;
		mysmoother = @HMMsmoother;
	end
	
	n = 10000;
	tps = zeros(nbReps, 3);
	for rep = 1:nbReps, rep
	% create the sample
			tic; 
			[x,y] = HMMsample(nu, Q, g, n);
			tps(rep, 1) = toc;

	% filtering / smoothing    
			tic; 
			[phi, c] = myfilter(y, nu, Q, g);
			beta = mysmoother(y, Q, g, c);
			post = phi .* beta;
			tps(rep, 2) = toc;

	% estimate parameters from observations only
			nbTrials = 10; 
			bestl=-inf;
			tic; 
			for j=1:nbTrials
				[Qn, gn, l] = HMMbaumwelch(y, nu); 
				if l>bestl
					Qh = Qn; gh = gn; bestl = l;
				end
			end
			tps(rep, 3) = toc;
	end
	tps
	mean(tps, 1)
	if finish, exit; end
