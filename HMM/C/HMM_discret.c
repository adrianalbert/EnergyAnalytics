/*
 *  written by Aurelien Garivier, CNRS & Telecom Paristech
 *  January 2012
 * (small bug corrected by Julius Su, Caltech)
 *
 * Baum-Welch algorithm for discrete Hidden Markov models
 * see http://www.telecom-paristech.fr/~garivier/code/index.html
 * 
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

/* scenario 1*//*
#define N 2
#define M 4
*/

/* scenario 2-3*/
#define N 50
#define M N
#define epsilon 0.1

	       
#define n 10000
#define nbReps 3


/*
sample from discrete distribution
in:   p = vector of probabilities,assumed to sum to 1
*/
int randm(double p[N]){
	int res=0;
	double q=p[0];
	double u=(rand()+0.0)/RAND_MAX;
	while(u>q) q+=p[++res];
	return(res);
}

/*
sample a trajectory from a hidden markov chain
in:   nu = initial distribution as vector of size k
      Q = transition matrix of size k
      n = positive integer
out:  (x,y) = sample trajectory of size n of a HMM defined by (nu, Q, g):
      x = sample trajectory of size n of a Markov Chain with initial distribution nu and transition matrix Q
      y = observations such that the conditionnal distribution of y[k]
      given x[k] is g(x[k], :)
*/
void HMMsample(double nu[N], double Q[N][N], double g[N][M], int x[n], int y[n]){
	int k;
	x[0] = randm(nu);
	y[0] = randm(g[x[0]]);
	for(k=1; k<n; k++){
		x[k] = randm(Q[x[k-1]]);
		y[k] = randm(g[x[k]]);
	}
} 

/*
HMM filtering of an observation sequence, given hmm parameters
in:   y = vector of observations, assumed to be in range(g.shape[1])
      nu = initial distribution as vector of size k
      Q = transition matrix of size k x k
      g = emission matrix with k rows
out:  phi = filter: P(x[t]=x | y[0:t]=y[0:t]) for 0<=x<k and 0<=t<n
      c(t) = conditional likelihood: P(Y[t] = y[t]| Y[0:t-1]=y[0:t-1])
*/
void HMMfilter(int y[n], double nu[N], double Q[N][N], double g[N][M], double phi[n][N], double c[n]){
	int i,j,t;
	double z[N];
	c[0]=0;
	for(j=0; j<N; j++){
		z[j] = nu[j]*g[j][y[0]];
		c[0] += z[j];
	}
	for(j=0; j<N; j++) phi[0][j] = z[j]/c[0];
	for(t=1; t<n; t++){
		c[t]=0;
		for(j=0; j<N; j++){
			z[j]=0;
			for(i=0; i<N; i++) z[j]+=phi[t-1][i]*Q[i][j]*g[j][y[t]];
			c[t] += z[j];
		}
		for(j=0; j<N; j++) phi[t][j] = z[j]/c[t];
	}
}

/*
HMM filtering of an observation sequence, given hmm parameters
in:   y = vector of observations, assumed to be in range(Q.shape[0])
      Q = transition matrix of size k x k
      g = emission matrix with k rows
      c = conditional likelihoods, computed by HMMfilter
out:  beta = smoothing factors: P(y[t+1:n]=y[t+1:n] | X[t]=x) / P(Y[t+1:n]=y[t+1:n] | Y[0:t]=y[1:t]) for 0<=x<k and 1<=t<n
      permits to compute the posterior distribution of the hidden states 
      P(X[t]=x | Y[0:n]=y[0:n])  as post = phi .* beta
*/
void HMMsmoother(int y[n], double Q[N][N], double g[N][M], double c[n], double beta[n][N]){
	int i,j,t;
	for(j=0; j<N; j++) beta[n-1][j] = 1;
	for(t=n-2; t>=0; t--){
		for(i=0; i<N; i++){
			double z=0;
			for(j=0; j<N; j++)
				z += Q[i][j] * g[j][y[t+1]] * beta[t+1][j];
			beta[t][i] = z / c[t+1];
		}
	}
}

/*
utility functions: sample random transition and emission kernels
*/
void randomTransitionKernel(double K[N][N]){
	int i,j;
	double s;
	for(i=0; i<N; i++){
		s=0;
		for(j=0; j<N; j++)
			s+=(K[i][j] = (rand()+0.0)/RAND_MAX);
		for(j=0; j<N; j++)
			K[i][j] /= s;
	}
}

void randomEmissionKernel(double K[N][M]){
	int i,j;
	double s;
	for(i=0; i<N; i++){
		s=0;
		for(j=0; j<M; j++)
			s+=(K[i][j] = (rand()+0.0)/RAND_MAX);
		for(j=0; j<M; j++)
			K[i][j] /= s;
	}
}

/*
compute maximum likehood estimate using Expectation-Maximization
iterations
in:   y = vector of observations 
      nu = initial distribution of the hidden chain
      tol = tolerance for the stopping criterion
      maxIt = maximal number of iterations
out:  Q = estimate of the transition matrix of the hidden markov process
      g = estimated probabilities of transition: g(x,y) = estimate of P(Y=y | X=x) for 0<=x<k
      l = log-likelihood of y for parameters Q and g
*/
double HMMbaumwelch(int y[n], double nu[N], double Q[N][N], double g[N][M]){
	// returns best log-likelihood found
	static const int maxIt = 100;
	static const double tol = 1e-4;
	int i, j, it, t;
	double z, l=0, change=tol+1;
	randomTransitionKernel(Q);
	randomEmissionKernel(g);
	static double phi[n][N];
	static double beta[n][N];
	static double c[n];
	double A[N][M];
	double s[N];
	double B[N][N];
	double s2[N];
	for(it=0; (change > tol) && (it<maxIt); it++){
		change = 0;
		HMMfilter(y, nu, Q, g, phi, c);
		HMMsmoother(y, Q, g, c, beta);
		for(i=0; i<N; i++){
			s[i]=0; s2[i]=0;
			for(j=0; j<M; j++) A[i][j] = 0;
			for(j=0; j<N; j++) B[i][j] = 0;
		}
		for(t=0; t<n; t++)
			for(i=0; i<N; i++){
				z=phi[t][i]*beta[t][i];
				A[i][y[t]] += z;
				s[i] += z;
			}
		for(i=0; i<N; i++)
			for(j=0; j<M; j++){
				change = MAX(change, fabs(g[i][j] - A[i][j] / s[i]));
				g[i][j] = A[i][j] / s[i];
			}
		for(t=1; t<n; t++){			
			for(i=0; i<N; i++)
				for(j=0; j<N; j++){
					z=phi[t-1][i]*Q[i][j]*g[j][y[t]]*beta[t][j]/c[t]; 
					B[i][j] += z;
					s2[i] += z;
				}
		}
		for(i=0; i<N; i++)
			for(j=0; j<N; j++){
				change = MAX(change, fabs(Q[i][j] - B[i][j] / s2[i]));
				Q[i][j] = B[i][j] / s2[i];
			}
	}
	for(t=0; t<n; t++) l += log(c[t]);
	return l;
}

/*
BENCHMARK for the comparison with other programming languages
http://www.telecom-paristech.fr/~garivier/code/index.html
*/
int main(){
	int i,j,rep,t, nbTrials = 10;
	double l, bestl;
	static int x[n]; 
	static int y[n]; 
	static double c[n];
	static double phi[n][N];
	static double beta[n][N];
	double Qh[N][N];
	double gh[N][M];
	clock_t begin, end;
	double diffms;
	double tps[nbReps][3];
	double avg[3];
	/* scenario 1: *//*
	double nu[N] = {0, 1};
	double Q[N][N] = {{0.8, 0.2},{0.1, 0.9}};
	double g[N][M] = {{0.25, 0.25, 0.25, 0.25},{0.1, 0.1, 0.4, 0.4}};
	*/
	/* scenario 2-3:*/
	double nu[N]; 
	double Q[N][N];
	double g[N][M];
	nu[0] = 1; for(i=1; i<N; i++) nu[i] = 0;
	for (i=0; i<N; i++) for(j=0; j<N; j++)
		if (i==j) Q[i][j] = 1-2*epsilon;
		else if ((i==0 && j==1) || (i==N-1 && j==N-2)) Q[i][j] = 2*epsilon;
		else if (abs(i-j)==1) Q[i][j] = epsilon;		
		else Q[i][j] = 0;		
	for (i=0; i<N; i++) for(j=0; j<M; j++) g[i][j] = Q[i][j];
	
	
	srand(time(NULL));
	for(rep=0; rep<nbReps; rep++){	printf("%d / %d\n", rep, nbReps);
		begin=clock();
		for(j=0; j<100; j++) HMMsample(nu, Q, g, x, y);
		end=clock();
		
		tps[rep][0] = (end-begin)/(CLOCKS_PER_SEC/1000.0)/100;

		begin=clock();
		for(j=0; j<1000; j++){
			HMMfilter(y, nu, Q, g, phi, c);
			HMMsmoother(y, Q, g, c, beta);
			for(t=0; t<n; t++) for(i=0; i<N; i++) phi[t][i]*=beta[t][i];
		}
		end=clock();
		tps[rep][1] = (end-begin)/(CLOCKS_PER_SEC/1000.0)/1000;

		bestl = -1e300;
		begin=clock();
		for(i=0; i<nbTrials; i++){	
			l = HMMbaumwelch(y, nu, Q, g);
			if (l>bestl){
				bestl = l;
				for(i=0; i<N; i++) for(j=0; j<N; j++) Qh[i][j] = Q[i][j];
				for(i=0; i<N; i++) for(j=0; j<M; j++) gh[i][j] = g[i][j];
			}
		}
		end=clock();
		tps[rep][2] = (end-begin)/(CLOCKS_PER_SEC / 1000.0);

	}
	
	printf("\n\ntimer:\n");
	for(rep=0; rep<nbReps; rep++){
		for (j=0; j<3; j++){
			printf("%.3f ", tps[rep][j]);
			avg[j] += tps[rep][j]/nbReps;
		}
		printf("\n");
	}
	printf("\naverage :\nsample: %.3f \nfiltering/smoothing: %.3f \nEM: %.3f\n", avg[0], avg[1], avg[2]);
	
	return(0);
}
