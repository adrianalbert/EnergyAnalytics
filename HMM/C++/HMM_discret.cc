/*
 * HMMarma.cpp
 *
 *  Created on: 3 févr. 2012
 *      Author: Pierre Pudlo
 *     (modified by: Aurelien Garivier)
 *
 *  Another way of programming the Baum-Welch algorithm for HMM
 *  with Armadillo, see http://arma.sourceforge.net/
 *
 *  Compile with:
 *  g++ -O3 HMMarma.cpp -o HMMarma -larmadillo
 *
 */

// disable run-time bound checks
#define ARMA_NO_DEBUG
#include <armadillo>
#include <limits>
#include <utility>
#include <cstdlib>
#include <ctime>

using namespace arma;


class HMMsample
{
public:
    void sample(vec & nu, mat & Q, mat & g, int n);
    ivec _x; // hidden states
    ivec _y; // observations
};

void HMMsample::sample(vec & nu, mat & Q, mat & g, int n){
	mat cQ = cumsum(Q, 1);
	mat cg = cumsum(g, 1);

	_x = zeros<ivec>(n);
	_y = zeros<ivec>(n);


	double random = (double)rand()/RAND_MAX;
	_x(0) = as_scalar(accu(cumsum(nu) < random ));
	random = (double)rand()/RAND_MAX;
	_y(0) = as_scalar(accu(cg.row(_x(0)) < random));
	for(size_t j=1; j<n; ++j){
		random = (double)rand()/RAND_MAX;
		_x(j) = as_scalar(accu(cQ.row(_x(j-1)) < random));
		random = (double)rand()/RAND_MAX;
		_y(j) = as_scalar(accu(cg.row(_x(j)) < random));
	}
}

class HMMfilter{
public:
	void compute(ivec &y, vec &nu, mat &Q, mat &g);
	mat _phi;
	vec _c;
};

void HMMfilter::compute(ivec &y, vec &nu, mat &Q, mat &g){
	int n = y.n_rows;
	_phi = zeros<mat>(Q.n_rows, n);
	_c = zeros<vec>(n);

	vec Z = nu % g.col(y(0));
	_c(0) = accu(Z);
	_phi.col(0) = Z / _c(0);

	for(size_t t = 1; t < n; ++t){
		Z = trans(trans(_phi.col(t-1)) * Q) % g.col(y(t));
		_c(t) = accu(Z);
		_phi.col(t) = Z / _c(t);
	}
}


class HMMsmoother{
public:
	void compute(ivec &y, mat &Q, mat &g, vec &c);
	mat _betaa;
};

void HMMsmoother::compute(ivec &y, mat &Q, mat &g, vec &c){
	int n = y.n_rows;
	_betaa = ones<mat>(Q.n_rows, n);
	for(int t=n-2; t>=0; --t){
		_betaa.col(t) = Q * (g.col(y(t+1)) % _betaa.col(t+1)) / c(t+1);
	}
}


class HMMbaumwelch{
public:
  void compute(ivec &y, vec &nu, double tol=1e-4, int maxIt = 100);
  mat _Q;
  mat _g;
  double _l;
};

void  HMMbaumwelch::compute(ivec &y, vec &nu, double tol, int maxIt){
  size_t k = nu.n_rows;
  size_t r = max(y)+1;
  size_t n = y.n_rows;
  imat Y = zeros<imat>(n, r);
  for(size_t i=0; i<n; ++i) Y(i, y(i)) = 1;

  _Q = randu<mat>(k, k);
  _Q = _Q / (sum(_Q, 1) * ones<mat>(1,k));
  _g = randu<mat>(k, r);
  _g = _g / (sum(_g,1) * ones<mat>(1,r));

  double it = 0;
  mat oldQ = _Q;
  mat oldg = _g + tol + 1;
  HMMfilter my_filter;
  HMMsmoother my_smoother;
  mat gaty = zeros<mat>(k, n-1);
  while( (norm(oldQ - _Q, 1) + norm(oldg - _g, 1) > tol) && (it < maxIt)){
	  ++it;
	  // compute the posterior distribution for the current parameters

	  my_filter.compute(y, nu, _Q, _g);

	  my_smoother.compute(y, _Q, _g, my_filter._c);
	  mat post = my_filter._phi % my_smoother._betaa;

	  // expectation of the number of transitions under the current parameters

	  for(size_t j=0; j<n-1; ++j)
		  gaty.col(j) = _g.col(y(j+1));

	  mat N = _Q % (
			  my_filter._phi.cols(0, n-2) * trans(
					  my_smoother._betaa.cols(1, n-1) % gaty / (ones<mat>(k,1) * trans(my_filter._c.subvec(1,n-1)))
					  )) ;
	  //cout << "N = "<< N << endl;
	  // expectation of the number of emissions
	  mat M = post * Y;
	  //cout << "M = " << M << endl;

	  // re-estimation
	  oldQ = _Q; oldg = _g;
	  _Q = N / (sum(N,1) * ones<mat>(1, k));
	  _g = M / (sum(M,1) * ones<mat>(1, r));
  }
  _l = accu(log(my_filter._c));
}



int main()
{
	srand(time(NULL));
	//scenario 1 
	
	vec nu = "0. 1.";
	mat Q = "0.8 0.2; 0.1 0.9";
	mat g = "0.25 0.25 0.25 0.25; 0.1 0.1 0.4 0.4";
	
	//scenario 2
	/*
	int N = 7;
	double epsilon=0.1;
	vec nu = zeros<vec>(N);
	nu(1)=1;
	mat Q = zeros<mat>(N, N);
	for (int i=0; i<N; i++) for(int j=0; j<N; j++)
		if (i==j) Q(i,j) = 1-2*epsilon;
		else if ((i==0 && j==1) || (i==N-1 && j==N-2)) Q(i,j) = 2*epsilon;
		else if (abs(i-j)==1) Q(i,j) = epsilon;		
		else Q(i,j) = 0;
	mat g=Q;
	*/
	int n = 10000;
	int nbReps = 10;
	mat tps = zeros<mat>(nbReps,3);
	clock_t begin, end;

	for(int rep = 0; rep < nbReps; ++rep){
		cout << "Replicate number " << rep+1 << " over " << nbReps << endl;
		HMMsample my_sample;
		begin = clock();
		for(int j = 0; j<100; ++j)
			my_sample.sample(nu, Q, g, n);
		end = clock();
		tps(rep, 0) = (end - begin) / (CLOCKS_PER_SEC / 1000.0) / 100;
		HMMfilter my_filter;
		HMMsmoother my_smoother;
		begin = clock();
		for(int j = 0; j < 1000; ++j){
			my_filter.compute(my_sample._y, nu, Q, g);
			my_smoother.compute(my_sample._y, Q, g, my_filter._c);
			mat post = my_filter._phi % my_smoother._betaa;
		}
		end = clock();
		tps(rep, 1) = (end - begin) / (CLOCKS_PER_SEC / 1000.0) / 1000;
		// double l0 = accu(log(my_filter._c)); // log-likelihood of the observation
		// cout << "Log-likelihood of the observation : " << l0 << endl;

		// estimate parameters from observations only
		int nbTrials = 10;
		double bestl = -1e30;
		begin = clock();
		HMMbaumwelch mybw;
		HMMbaumwelch best;
		for(int j = 0; j < nbTrials; ++j){
			mybw.compute(my_sample._y, nu);
			if( mybw._l > bestl) {
				best = mybw;
				bestl = mybw._l;
			}
		}
		end = clock();
		tps(rep, 2) = (end - begin) / (CLOCKS_PER_SEC / 1000.0);
		//cout << "Likelihood = " << bestl << endl;
		//cout << "hat Q = " << best._Q << endl;
		//cout << "hat g =" << best._g << endl;
	}

	cout << "Timer:" << endl;
	cout << tps << endl;

	cout << "Average :" << endl;
	cout << mean(tps) << endl;

	return 0;
}
