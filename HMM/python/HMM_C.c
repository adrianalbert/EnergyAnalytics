/*
 *  written by Aurelien Garivier, CNRS & Telecom Paristech
 *  January 2012
 *
 * Baum-Welch algorithm for discrete Hidden Markov models
 * see http://www.telecom-paristech.fr/~garivier/code/index.html
 * 
 * *=================================================================
 *
 * compile with "make", or:
 * mex -v  HMMfilter.c
 * optimized: mex COPTIMFLAGS='-O2' -v HMMfilter.c 
 *
 */


#include <Python.h>
#include <math.h>
#include <numpy/arrayobject.h>

void _HMMfilter(int N, int M, int n, double* y, double* nu, double* Q, double* g, double* phi, double* c){
	int i,j,t;
	double z[N];
	c[0]=0;
	for(j=0; j<N; j++){
	  z[j] = nu[j]*g[j*M+((int)y[0])];
		c[0] += z[j];
	}
	for(j=0; j<N; j++) phi[j*n+0] = z[j]/c[0];
	for(t=1; t<n; t++){
		c[t]=0;
		for(j=0; j<N; j++){
			z[j]=0;
			for(i=0; i<N; i++) z[j]+=phi[i*n+(t-1)]*Q[i*N+j]*g[j*M+((int)y[t])];
			c[t] += z[j];
		}
		for(j=0; j<N; j++) phi[j*n + t] = z[j]/c[t];
	}
}

/*  (phi, c) = HMMfilter(y, nu, Q, g):
in:   y = vector of observations of size n, with values between 1 and r
      nu = initial distribution as vector of size k
      Q = transition matrix of size k
      g = emission matrix of size k x r
out:  phi = filter(x,t) = P(X(t)=x | Y(1:t)=y(1:t)) for 1<=x<=k and 1<=t<=n
c(t) = P(Y(t) = y(t)| Y(1:t-1)=y(1:t-1))*/

static PyObject* HMMfilter(PyObject* self, PyObject* args)
{
  int N,M,n;
  double *y, *nu, *Q, *g, *phi, *c;
  PyObject *PY_y, *PY_nu, *PY_Q, *PY_g, *PY_phi, *PY_c;
  PyObject *PYs_y, *PYs_nu, *PYs_Q, *PYs_g, *PY_out;
  npy_intp phiSize[2];

  if (!PyArg_ParseTuple(args, "OOOO", &PY_y, &PY_nu, &PY_Q, &PY_g))
    return NULL;
        
  /* type and contigency checking */
  PYs_y = PyArray_FROM_OTF(PY_y, NPY_DOUBLE, NPY_IN_ARRAY);
  PYs_nu = PyArray_FROM_OTF(PY_nu, NPY_DOUBLE, NPY_IN_ARRAY);
  PYs_Q = PyArray_FROM_OTF(PY_Q, NPY_DOUBLE, NPY_IN_ARRAY);
  PYs_g = PyArray_FROM_OTF(PY_g, NPY_DOUBLE, NPY_IN_ARRAY);
    
  N = PyArray_DIMS(PYs_Q)[0];    
  M = PyArray_DIMS(PYs_g)[1];
  n = PyArray_DIMS(PYs_y)[0];
  y = (double *)PyArray_DATA(PYs_y);
  nu = (double *)PyArray_DATA(PYs_nu);
  Q = (double *)PyArray_DATA(PYs_Q);
  g = (double *)PyArray_DATA(PYs_g);

  phiSize[0] = N; phiSize[1] = n;
  PY_phi = PyArray_SimpleNew(2, phiSize, PyArray_DOUBLE);
  PY_c = PyArray_SimpleNew(1, &phiSize[1], PyArray_DOUBLE);
  phi = (double *)PyArray_DATA(PY_phi);
  c = (double*)PyArray_DATA(PY_c);
  PY_out = Py_BuildValue("(O,O)", PY_phi, PY_c);

  _HMMfilter(N,M,n,y,nu,Q,g,phi,c);
    
  Py_DECREF(PYs_y); Py_DECREF(PYs_nu); Py_DECREF(PYs_Q); Py_DECREF(PYs_g);
  Py_DECREF(PY_phi); Py_DECREF(PY_c);
  return PY_out; 
}



void _HMMsmoother(int N, int M, int n, double* y, double* Q, double* g, double* c, double* beta){
	int i,j,t;
	for(j=0; j<N; j++) beta[j*n+(n-1)] = 1;
	for(t=n-2; t>=0; t--){
		for(i=0; i<N; i++){
			double z=0;
			for(j=0; j<N; j++)
			  z += Q[i*N+j] * g[j*M+((int)y[t+1])] * beta[j*n+(t+1)];
			beta[i*n+t] = z / c[t+1];
		}
	}
}


/* beta = smoother(y, Q, g, c): */
static PyObject* HMMsmoother(PyObject* self, PyObject* args)
{
  int N,M,n;
  double *y, *Q, *g, *c, *beta;
  PyObject *PY_y, *PY_Q, *PY_g, *PY_c, *PY_beta;
  PyObject *PYs_y, *PYs_Q, *PYs_g, *PYs_c;
  npy_intp phiSize[2];

  if (!PyArg_ParseTuple(args, "OOOO", &PY_y, &PY_Q, &PY_g, &PY_c))
    return NULL;
        
  /* type and contigency checking */
  PYs_y = PyArray_FROM_OTF(PY_y, NPY_DOUBLE, NPY_IN_ARRAY);
  PYs_Q = PyArray_FROM_OTF(PY_Q, NPY_DOUBLE, NPY_IN_ARRAY);
  PYs_g = PyArray_FROM_OTF(PY_g, NPY_DOUBLE, NPY_IN_ARRAY);
  PYs_c = PyArray_FROM_OTF(PY_c, NPY_DOUBLE, NPY_IN_ARRAY);
    
  N = PyArray_DIMS(PYs_Q)[0];    
  M = PyArray_DIMS(PYs_g)[1];
  n = PyArray_DIMS(PYs_y)[0];
  y = (double *)PyArray_DATA(PYs_y);
  Q = (double *)PyArray_DATA(PYs_Q);
  g = (double *)PyArray_DATA(PYs_g);
  c = (double *)PyArray_DATA(PYs_c);

  phiSize[0] = N; phiSize[1] = n;
  PY_beta = PyArray_SimpleNew(2, phiSize, PyArray_DOUBLE);
  beta = (double *)PyArray_DATA(PY_beta);
  c = (double*)PyArray_DATA(PY_c);

  _HMMsmoother(N,M,n,y,Q,g,c,beta);
    
  Py_DECREF(PYs_y); Py_DECREF(PYs_Q); Py_DECREF(PYs_g); Py_DECREF(PYs_c);

  return PY_beta; 
}

 
static PyMethodDef HMM_CMethods[] = {
  {"HMMfilter", HMMfilter, METH_VARARGS, "HMMfilter"},
  {"HMMsmoother", HMMsmoother, METH_VARARGS, "HMMsmoother"},
  {NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initHMM_C(void)
{
  (void) Py_InitModule("HMM_C", HMM_CMethods);
  import_array();
}
