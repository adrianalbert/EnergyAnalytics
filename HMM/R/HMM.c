#include <R.h>
#include <stdio.h>
#include <math.h>


void HMMfilter(int* N, int* M, int* n, double* y, double* nu, double* Q, double* g, double* phi, double* c){
	int i,j,t;
	double z[*N];
	c[0]=0;
	for(j=0; j<*N; j++){
		 z[j] = nu[j]*g[j+(*N)*((int)y[0]-1)];
		 c[0] += z[j];
	}
	for(j=0; j<*N; j++) phi[j+0*(*N)] = z[j]/c[0];
	for(t=1; t<*n; t++){
		c[t]=0;
		for(j=0; j<*N; j++){
			z[j]=0;
			for(i=0; i<*N; i++) z[j]+=phi[i+(*N)*(t-1)]*Q[i+(*N)*j]*g[j+(*N)*((int)y[t]-1)];
			c[t] += z[j];
		}
		for(j=0; j<(*N); j++) phi[j + *(N)*t] = z[j]/c[t];
	}
}


void HMMsmoother(int* N, int* M, int* n, double* y, double* Q, double* g, double* c, double* beta){
	int i,j,t;
	for(j=0; j<*N; j++) beta[j+(*N)*((*n)-1)] = 1;
	for(t=(*n)-2; t>=0; t--){
		for(i=0; i<*N; i++){
			double z=0;
			for(j=0; j<*N; j++)
			  z += Q[i+(*N)*j] * g[j+(*N)*((int)y[t+1] -1)] * beta[j+(*N)*(t+1)];
			beta[i+(*N)*t] = z / c[t+1];
		}
	}
}
