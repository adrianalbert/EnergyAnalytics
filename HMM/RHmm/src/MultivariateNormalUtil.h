/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: MultivariateNormalUtil.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _MULTIVARIATENORMALUTIL_H_
#define _MULTIVARIATENORMALUTIL_H_
#pragma once
#include "OTMathUtil.h"

#ifndef _RDLL_
        #include "REquivalents.h"
#else
        #include <cstdio>
        #include <cstdlib>
        #include <R.h>
        #include <Rinternals.h>
        #include <Rmath.h>
#endif // _RDLL_


#ifndef SQRT_TWO_PI
        #define SQRT_TWO_PI     2.5066282746310002
#endif //SQRT_TWO_PI
#ifndef uint
        typedef unsigned int uint ;
#endif //int

extern void SymetricInverseAndDet(cDMatrix& theMat, double& theDet, cDMatrix& theInvMat) ;

extern void MultivariateNormalDensity(cDVector& thex, cDVector& theMu, cDMatrix& theInvCov, double theDet, double* theDens) ;

extern void MultivariateNormalDensity(cDVector& thex, cDVector& theMu, cDMatrix& theInvCov, double theDet, cDVector& theDens) ;

extern void SymDetDeriv(cDMatrix& theMat, cDVector& theGrad, cDMatrix& theHess) ;

extern void InvCovMatDeriv(cDMatrix& theInvCov, cDMatrix* theGrad, cDMatrix** theHess) ;

extern void MultivariateNormalDensityDeriv(cDVector& thex, cDVector& theMu, cDMatrix& theCov, cDMatrix& theInvCov, double theDet, cDVector& theDensity, cDVector* theGrad, cDMatrix* theHess) ;

extern void MultivariateNormalDensityDeriv(cDVector& thex, cDVector& theMu, cDMatrix& theCov, cDMatrix& theInvCov, double theDet, cDVector* theGrad, cDMatrix* theHess) ;

#endif //_MULTIVARIATENORMALUTIL_H_
