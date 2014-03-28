/**************************************************************
 *** RHmm package
 ***                                                         
 *** File: LogProb.h 
 ***                                                         
 *** Author: Ollivier TARAMASCO <Ollivier.Taramasco@imag.fr> 
 *** Author: Sebastian BAUER <sebastian.bauer@charite.de>
 ***                                                         
 **************************************************************/

#ifndef _LOGPROB_H_
#define _LOGPROB_H_
#include <cfloat>
#include <cmath>
#include "math.h"
        #define LOGZERO -DBL_MAX
#ifndef MIN
        #define MIN(p,q) ((p) < (q) ? (p) : (q))
#endif //MIN
#ifndef MAX
        #define MAX(p,q) ((p) > (q) ? (p) : (q))
#endif //MAX
extern double eexp(const double theX) ;
extern double eln(const double theX) ;
extern double elnsum(const double theX, const double theY) ;
extern double elnproduct(const double theX, const double theY) ;

#endif // _LOGPROB_H_


